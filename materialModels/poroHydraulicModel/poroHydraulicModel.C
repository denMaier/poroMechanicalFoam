/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "poroHydraulicModel.H"
#include "volFields.H"


namespace Foam{
    defineTypeNameAndDebug(poroHydraulicModel, 0);
}

namespace
{
    Foam::word lookupLawName
    (
        const Foam::dictionary& zoneSubDict,
        const Foam::word& canonicalKey,
        const Foam::word& legacyKey,
        const Foam::word& zoneName
    )
    {
        if (zoneSubDict.found(canonicalKey))
        {
            if (zoneSubDict.found(legacyKey))
            {
                WarningInFunction
                    << "Zone '" << zoneName << "' defines both '" << canonicalKey
                    << "' and legacy key '" << legacyKey << "'. Using '"
                    << canonicalKey << "'."
                    << Foam::endl;
            }

            return zoneSubDict.get<Foam::word>(canonicalKey);
        }

        if (zoneSubDict.found(legacyKey))
        {
            WarningInFunction
                << "Zone '" << zoneName << "' uses legacy key '" << legacyKey
                << "'. Please rename it to '" << canonicalKey
                << "' in poroHydraulicProperties."
                << Foam::endl;

            return zoneSubDict.get<Foam::word>(legacyKey);
        }

        FatalIOErrorInFunction(zoneSubDict)
            << "Zone '" << zoneName << "' does not define required selector key '"
            << canonicalKey << "'"
            << " (legacy alias '" << legacyKey << "' is also accepted)."
            << Foam::exit(Foam::FatalIOError);

        return Foam::word::null;
    }

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
}

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
poroHydraulicModel::poroHydraulicModel
(
    const volScalarField& pField,
    const dimensionedVector& gravity
)
    : regIOobject
    (
        IOobject
        (
            "poroHydraulicModel",
            pField.time().constant(),
            pField.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    poroHydraulicProperties_
    (
        IOobject
        (
            "poroHydraulicProperties",
            pField.time().constant(),
            pField.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    pField_(pField),
    storageLawPtr_(pField.mesh().cellZones().size()),
    conductivityModelPtr_(pField.mesh().cellZones().size()),
    rho_
    (
        poroHydraulicProperties_.lookupOrAddDefault<dimensionedScalar>
        (
            "rho",
            dimensionedScalar("rho", dimDensity, 1000)
        )
    ),
    gamma_
    (
        IOobject
        (
            "gamma_water",
            mesh().time().constant(),
            pField.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*gravity
    ),
    href_
    (
        poroHydraulicProperties_.lookupOrAddDefault<dimensionedScalar>
        (
            "href",
            dimensionedScalar("href", dimLength, 0.0)
        )
    ),
    z_
    ( //  geodetic height (depends on direction of gravity)
        IOobject
        (
            "z",
            mesh().time().constant(),
            pField.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        - mesh().C() & vector(gamma_.value()).normalise()
    ),
    p_Hyd_
    ( //  hydrostatic pressure
        IOobject
        (
            "p_Hyd",
            mesh().time().timeName(),
            pField.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (href_ - z_) * magGamma()
    ),
    updatesSs_(false),
    updatesK_(false)
{

    Info << "Creating the poroHydraulicModel" << nl
         << "Gravity direction: " << vector(gamma_.value()).normalise() << nl
         << "Water specific weight: " << magGamma() << nl
	 << "Referential Watertable is at z = " << href_ << endl;

    // Build one constitutive law pair per cellZone. The model assumes that
    // every hydraulic region is represented by a cellZone with its own sub-dict.
    forAll(pField.mesh().cellZones(),zoneI)
    {
        const word& zoneName = pField.mesh().cellZones()[zoneI].name();
        const dictionary& zoneSubDict =
             poroHydraulicProperties_.optionalSubDict 
                      (
                        zoneName
                      );

        storageLawPtr_.set
        (
            zoneI,
            storageLaw::New
            (
                lookupLawName(zoneSubDict, "storageLaw", "StorageModel", zoneName),
                // Storage-law constructors may insert default entries so a
                // later `_withDefaultValues` write-back reflects the effective setup.
                const_cast<dictionary&>(zoneSubDict),
                pField_
            )
        );

        updatesSs_ = (updatesSs_ || storageLawPtr_[zoneI].updatesSs());

        conductivityModelPtr_.set
        (
            zoneI,
            conductivityModel::New
            (
                zoneSubDict.found("conductivityLaw")
              ? lookupLawName(zoneSubDict, "conductivityLaw", "conductivityModel", zoneName)
              : zoneSubDict.found("conductivityModel")
                    ? lookupLawName(zoneSubDict, "conductivityLaw", "conductivityModel", zoneName)
                    : word("constant"),
                // Conductivity-model constructors may insert default entries so
                // a later `_withDefaultValues` write-back reflects the effective setup.
                const_cast<dictionary&>(zoneSubDict),
                pField
            )
        );

        updatesK_ = (updatesK_ || conductivityModelPtr_[zoneI].updatesK());
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

poroHydraulicModel::~poroHydraulicModel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> poroHydraulicModel::n0() const
{
    IOobject nHeader
    (
        "n",
        mesh().time().constant(),
        pField().db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    // Prefer a fully specified porosity field from disk when available.
    if (nHeader.typeHeaderOk<volScalarField>(true))
    {
        tmp<volScalarField> tn0
        (
            new volScalarField
            (
                nHeader,
                mesh()
            )
        );

        return tn0;
    }
    else
    {
        // Otherwise build a piecewise-constant porosity field from the
        // per-cellZone `n` entries in poroHydraulicProperties.
        tmp<volScalarField> tn0
        (
            new volScalarField
            (
                IOobject
                (
                    "n",
                    mesh().time().constant(),
                    pField().db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimless,0.0)
            )
        );

        volScalarField& n = tn0.ref();
        
        const cellZoneMesh& cellZones = hydraulicZones();
        PtrList<dimensionedScalar> nList(cellZones.size());

        forAll(cellZones,iZone)
        {
            const word& zoneName = cellZones[iZone].name();
            const dictionary& zoneSubDict =
                    poroHydraulicProperties().optionalSubDict 
                    (
                        zoneName
                    );
            nList.set
            (
                iZone,
                new dimensionedScalar(zoneSubDict.get<dimensionedScalar>("n"))
            );
        }

        assembleInternalScalarField
        (
            n,
            [&](const label zoneI, const label)
            {
                return nList[zoneI].value();
            }
        );

        assembleBoundaryScalarField
        (
            n,
            [&](const label zoneI, const label, const label)
            {
                return nList[zoneI].value();
            }
        );

        return tn0;
    }
}

tmp<volScalarField> poroHydraulicModel::Ss(const volScalarField& n, const volScalarField& p)
{
    tmp<volScalarField> tSs
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("Ss"),
            mesh(),
            dimensionedScalar(dimless / pField().dimensions() , 0.0),
            // Ss is a volumetric storage quantity. Boundary-face values are not
            // used as constitutive state, so a derived zeroGradient patch field
            // is sufficient here.
            "zeroGradient"
        )
    );

    volScalarField& Ss_ = tSs.ref();

    // Storage is evaluated cell-wise because Ss is only physically interpreted
    // as a per-volume quantity. Boundary values are not used to drive boundary
    // constitutive state, so we only assemble the internal field explicitly.
    assembleInternalScalarField
    (
        Ss_,
        [&](const label zoneI, const label cellI)
        {
            return storageLawPtr_[zoneI].Ss
            (
                n.internalField()[cellI],
                p.internalField()[cellI],
                cellI
            );
        }
    );

    Ss_.correctBoundaryConditions();

    return tSs;
}

const tmp<volScalarField> poroHydraulicModel::k() const
{
    tmp<volScalarField> tk
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("k"),
            mesh(),
            dimensionedScalar(dimLength / dimTime, 0.0)
        )
    );

    volScalarField& k_ = tk.ref();

    // Conductivity is evaluated per cell and then copied to boundary faces
    // from the owner-cell zone to keep zone-wise material behavior consistent.
    assembleInternalScalarField
    (
        k_,
        [&](const label zoneI, const label cellI)
        {
            return conductivityModelPtr_[zoneI].k(cellI);
        }
    );

    assembleBoundaryScalarField
    (
        k_,
        [&](const label zoneI, const label patchI, const label faceI)
        {
            return conductivityModelPtr_[zoneI].k(patchI, faceI);
        }
    );
    
    return tk;
}

const tmp<surfaceScalarField> poroHydraulicModel::kf() const
{
    tmp<volScalarField> tk = this->k();
    return fvc::interpolate(tk());
}

} // End of namespace Foam
// ************************************************************************* //
