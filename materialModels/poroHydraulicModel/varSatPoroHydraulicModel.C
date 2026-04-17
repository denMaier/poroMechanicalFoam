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

#include "varSatPoroHydraulicModel.H"
#include "volFields.H"


namespace Foam{

    defineTypeNameAndDebug(varSatPoroHydraulicModel, 0);
}

namespace
{
    Foam::word lookupSaturationLawName
    (
        const Foam::dictionary& zoneSubDict,
        const Foam::word& zoneName
    )
    {
        if (zoneSubDict.found("saturationLaw"))
        {
            if (zoneSubDict.found("SWCC"))
            {
                Foam::WarningInFunction
                    << "Zone '" << zoneName << "' defines both 'saturationLaw'"
                    << " and legacy key 'SWCC'. Using 'saturationLaw'."
                    << Foam::endl;
            }

            return zoneSubDict.get<Foam::word>("saturationLaw");
        }

        if (zoneSubDict.found("SWCC"))
        {
            Foam::WarningInFunction
                << "Zone '" << zoneName << "' uses legacy key 'SWCC'."
                << " Please rename it to 'saturationLaw' in poroHydraulicProperties."
                << Foam::endl;

            return zoneSubDict.get<Foam::word>("SWCC");
        }

        Foam::FatalIOErrorInFunction(zoneSubDict)
            << "Zone '" << zoneName << "' does not define required selector key"
            << " 'saturationLaw' (legacy alias 'SWCC' is also accepted)."
            << Foam::exit(Foam::FatalIOError);

        return Foam::word::null;
    }

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
}

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
varSatPoroHydraulicModel::varSatPoroHydraulicModel(
    const volScalarField &pField,
    const dimensionedVector& gravity)
    : poroHydraulicModel(pField,gravity),
      saturationLawPtr_(pField.mesh().cellZones().size())
{
    // Build one saturation law per cellZone so Richards-type quantities can
    // vary by material region in the same way as storage and conductivity.
    forAll(pField.mesh().cellZones(),iZone)
    {
        const word& zoneName = pField.mesh().cellZones().names()[iZone];
        const dictionary& zoneSubDict =
             poroHydraulicProperties().optionalSubDict
                      (
                        zoneName
                      );

        saturationLawPtr_.set
        (
            iZone,
            saturationLaw::New
            (
                lookupSaturationLawName(zoneSubDict, zoneName),
                // Saturation-law constructors may insert default entries so a
                // later `_withDefaultValues` write-back reflects the effective setup.
                const_cast<dictionary&>(zoneSubDict),
                pField
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

varSatPoroHydraulicModel::~varSatPoroHydraulicModel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> varSatPoroHydraulicModel::pStar() const
{
    tmp<volScalarField> tpStar
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("pStar"),
            mesh(),
            dimensionedScalar(pField().dimensions(), 0.0),
            "zeroGradient" // Boundaryfield not important for pStar
        )
    );

    volScalarField& pStar_ = tpStar.ref();

    // pStar is assembled cell-wise because each zone may use a different SWCC.
    assembleInternalScalarField
    (
        pStar_,
        [&](const label zoneI, const label cellI)
        {
            return saturationLawPtr_[zoneI].pStar(cellI);
        }
    );

    pStar_.correctBoundaryConditions();
    
    return tpStar;
}

const tmp<volScalarField> varSatPoroHydraulicModel::kEff(const volScalarField &p)
{
    return k() * kr(p);
}

const tmp<surfaceScalarField> varSatPoroHydraulicModel::kEfff(const volScalarField &p)
{
    // The variably saturated face mobility is the saturated face conductivity
    // scaled by the interpolated relative conductivity.
    return kf()*fvc::interpolate(this->kr(p));
}

tmp<volScalarField> varSatPoroHydraulicModel::kr(const volScalarField &p)
{
    tmp<volScalarField> tkr
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("kr"),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    volScalarField& kr_ = tkr.ref();

    // kr is evaluated on cells and boundary faces because boundary conditions
    // may depend on the local relative conductivity state at the boundary.
    assembleInternalScalarField
    (
        kr_,
        [&](const label zoneI, const label cellI)
        {
            return saturationLawPtr_[zoneI].kr(p.internalField()[cellI], cellI);
        }
    );

    assembleBoundaryScalarField
    (
        kr_,
        [&](const label zoneI, const label patchI, const label faceI)
        {
            return saturationLawPtr_[zoneI].kr
            (
                p.boundaryField()[patchI][faceI],
                patchI,
                faceI
            );
        }
    );
    
    return tkr;
}

tmp<volScalarField> varSatPoroHydraulicModel::S(const volScalarField &p)
{
    tmp<volScalarField> tS
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("S"),
            mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    volScalarField& S_ = tS.ref();

    // Saturation follows the same per-zone evaluation pattern as kr because
    // boundary conditions may need the actual boundary saturation state.
    assembleInternalScalarField
    (
        S_,
        [&](const label zoneI, const label cellI)
        {
            return saturationLawPtr_[zoneI].S(p.internalField()[cellI], cellI);
        }
    );

    assembleBoundaryScalarField
    (
        S_,
        [&](const label zoneI, const label patchI, const label faceI)
        {
            return saturationLawPtr_[zoneI].S
            (
                p.boundaryField()[patchI][faceI],
                patchI,
                faceI
            );
        }
    );
    
    return tS;
}

tmp<volScalarField> varSatPoroHydraulicModel::C(const volScalarField &p)
{
    tmp<volScalarField> tC
    (
        new volScalarField
        (
            assembledScratchFieldIOobject("C"),
            mesh(),
            dimensionedScalar(dimless/pField().dimensions(), 0.0),
            // C is a volumetric dS/dp quantity used as cell-wise storage-like
            // information. Boundary-face values are not treated as constitutive
            // boundary state, so a derived zeroGradient patch field is enough.
            "zeroGradient"
        )
    );

    volScalarField& C_ = tC.ref();

    // C is only assembled on cells because it is interpreted per control
    // volume; nothing in the hydraulic model consumes an explicitly evaluated
    // boundary C state.
    assembleInternalScalarField
    (
        C_,
        [&](const label zoneI, const label cellI)
        {
            return saturationLawPtr_[zoneI].C(p.internalField()[cellI], cellI);
        }
    );

    C_.correctBoundaryConditions();
    
    return tC;
}



} // End of namespace Foam
// ************************************************************************* //
