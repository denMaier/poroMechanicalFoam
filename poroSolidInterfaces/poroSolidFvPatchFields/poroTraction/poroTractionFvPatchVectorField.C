/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "poroTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "poroFluidModel.H"
#include "mechanicalModel.H"
#include "poroMechanicalLaw2.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    poroTractionFvPatchVectorField::
        poroTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF)
        : solidTractionFvPatchVectorField(p, iF),
          totalTraction_(),
          tractionSeries_(),
          pressureSeries_(),
          secondOrder_(false),
          limitCoeff_(1.0),
          relaxFac_(1.0),
          pressureFieldName_("auto")
    {
        fvPatchVectorField::operator=(patchInternalField());
        gradient() = vector::zero;
    }

    poroTractionFvPatchVectorField::
        poroTractionFvPatchVectorField(
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const dictionary &dict)
        : solidTractionFvPatchVectorField(p, iF),
          totalTraction_(dict.lookupOrDefault<Switch>("total", true)),
          tractionSeries_(PatchFunction1<vector>::New(p.patch(), "traction", dict)),
          pressureSeries_(PatchFunction1<scalar>::New(p.patch(), "pressure", dict)),
          secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
          limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
          relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
          pressureFieldName_(dict.lookupOrDefault<word>("pressureFieldName", "auto"))
    {
        Info << "Creating " << type() << " boundary condition at patch " << patch().name() <<  nl;
        totalTraction_
        ? Info << tab << "Traction is total" << endl
        : Info << tab << "Traction is effective" << endl;

        if (!totalTraction_)
        {
            if
            (
                pressureFieldName_ != "auto"
             && pressureFieldName_ != "p"
             && pressureFieldName_ != "p_rgh"
            )
            {
                FatalIOErrorInFunction(dict)
                    << "Invalid pressureFieldName '" << pressureFieldName_
                    << "' on patch '" << patch().name() << "'." << nl
                    << "Valid values are 'auto', 'p', and 'p_rgh'."
                    << exit(FatalIOError);
            }

            if (dict.found("buoyancyIncluded"))
            {
                const Switch buoyancyIncluded
                (
                    dict.lookupOrDefault<Switch>("buoyancyIncluded", false)
                );

                const word legacyPressureName
                (
                    buoyancyIncluded ? "p_rgh" : "p"
                );

                if
                (
                    pressureFieldName_ != "auto"
                 && pressureFieldName_ != legacyPressureName
                )
                {
                    FatalIOErrorInFunction(dict)
                        << "Conflicting pressure selections on patch '"
                        << patch().name() << "': pressureFieldName is '"
                        << pressureFieldName_ << "', while buoyancyIncluded "
                        << buoyancyIncluded << " selects '"
                        << legacyPressureName << "'."
                        << exit(FatalIOError);
                }

                pressureFieldName_ = legacyPressureName;

                WarningInFunction
                    << "Entry 'buoyancyIncluded' on patch '" << patch().name()
                    << "' is deprecated. Use 'pressureFieldName auto' and let "
                    << "poroTraction inherit the pressure definition from "
                    << "poroMechanicalLaw2." << endl;
            }
        }

        if (dict.found("gradient"))
        {
            gradient() = vectorField("gradient", dict, p.size());
        }
        else
        {
            gradient() = vector::zero;
        }

        if (dict.found("value"))
        {
            Field<vector>::operator=(vectorField("value", dict, p.size()));
        }
        else
        {
            fvPatchVectorField::operator=(patchInternalField());
        }

        if (secondOrder_)
        {
            Info << "    second order correction" << endl;
        }

        if (limitCoeff_)
        {
            Info << "    limiter coefficient: " << limitCoeff_ << endl;
        }

        if (relaxFac_ < 1.0)
        {
            Info << "    relaxation factor: " << relaxFac_ << endl;
        }
    }

    poroTractionFvPatchVectorField::
        poroTractionFvPatchVectorField(
            const poroTractionFvPatchVectorField &stpvf,
            const fvPatch &p,
            const DimensionedField<vector, volMesh> &iF,
            const fvPatchFieldMapper &mapper)
        : solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
          totalTraction_(stpvf.totalTraction_),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          pressureFieldName_(stpvf.pressureFieldName_)
    {
    }

    poroTractionFvPatchVectorField::
        poroTractionFvPatchVectorField(
            const poroTractionFvPatchVectorField &stpvf)
        : solidTractionFvPatchVectorField(stpvf),
          totalTraction_(stpvf.totalTraction_),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          pressureFieldName_(stpvf.pressureFieldName_)
    {
    }

    poroTractionFvPatchVectorField::
        poroTractionFvPatchVectorField(
            const poroTractionFvPatchVectorField &stpvf,
            const DimensionedField<vector, volMesh> &iF)
        : solidTractionFvPatchVectorField(stpvf, iF),
          totalTraction_(stpvf.totalTraction_),
          tractionSeries_(stpvf.tractionSeries_.clone(this->patch().patch())),
          pressureSeries_(stpvf.pressureSeries_.clone(this->patch().patch())),
          secondOrder_(stpvf.secondOrder_),
          limitCoeff_(stpvf.limitCoeff_),
          relaxFac_(stpvf.relaxFac_),
          pressureFieldName_(stpvf.pressureFieldName_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void poroTractionFvPatchVectorField::autoMap(
        const fvPatchFieldMapper &m)
    {
        solidTractionFvPatchVectorField::autoMap(m);
    }

    // Reverse-map the given fvPatchField onto this fvPatchField
    void poroTractionFvPatchVectorField::rmap(
        const fvPatchVectorField &ptf,
        const labelList &addr)
    {
        solidTractionFvPatchVectorField::rmap(ptf, addr);
    }

    word poroTractionFvPatchVectorField::mechanicalPressureFieldName
    (
        const solidModel& solMod
    ) const
    {
        const mechanicalModel& mechanical = solMod.mechanical();
        const PtrList<mechanicalLaw>& laws = mechanical;

        word lawPressureName;
        bool foundPoroLaw = false;

        forAll(laws, lawI)
        {
            const poroMechanicalLaw2* poroLaw =
                dynamic_cast<const poroMechanicalLaw2*>(&laws[lawI]);

            if (!poroLaw)
            {
                continue;
            }

            if (!foundPoroLaw)
            {
                lawPressureName = poroLaw->pressureFieldName();
                foundPoroLaw = true;
            }
            else if (lawPressureName != poroLaw->pressureFieldName())
            {
                FatalErrorInFunction
                    << "poroTraction patch '" << patch().name()
                    << "' is attached to mechanical laws using different "
                    << "pressure fields ('" << lawPressureName << "' and '"
                    << poroLaw->pressureFieldName() << "'). A single traction "
                    << "patch cannot synchronize with both definitions."
                    << exit(FatalError);
            }
        }

        if (!foundPoroLaw)
        {
            if (pressureFieldName_ == "auto")
            {
                FatalErrorInFunction
                    << "Cannot select a pressure field automatically for "
                    << "poroTraction patch '" << patch().name() << "' because "
                    << "no poroMechanicalLaw2 was found. Set pressureFieldName "
                    << "explicitly to 'p' or 'p_rgh', or use "
                    << "poroMechanicalLaw2."
                    << exit(FatalError);
            }

            return pressureFieldName_;
        }

        if
        (
            pressureFieldName_ != "auto"
         && pressureFieldName_ != lawPressureName
        )
        {
            FatalErrorInFunction
                << "Pressure-field mismatch on poroTraction patch '"
                << patch().name() << "': the boundary condition selects '"
                << pressureFieldName_ << "', but poroMechanicalLaw2 uses '"
                << lawPressureName << "'. Use pressureFieldName auto or make "
                << "the two definitions consistent."
                << exit(FatalError);
        }

        return lawPressureName;
    }

    scalarField poroTractionFvPatchVectorField::patchBiotCoeff
    (
        const solidModel& solMod
    ) const
    {
        const mechanicalModel& mechanical = solMod.mechanical();
        const PtrList<mechanicalLaw>& laws = mechanical;

        if (laws.size() == 1)
        {
            const poroMechanicalLaw2* poroLaw =
                dynamic_cast<const poroMechanicalLaw2*>(&laws[0]);

            if (!poroLaw)
            {
                return scalarField(patch().size(), 1.0);
            }

            return scalarField
            (
                patch().size(),
                poroLaw->biotCoeffValue().value()
            );
        }

        tmp<volScalarField> tb
        (
            new volScalarField
            (
                IOobject
                (
                    "poroTractionBiotCoeff",
                    this->db().time().timeName(),
                    mechanical.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mechanical.mesh(),
                dimensionedScalar("one", dimless, 1.0),
                calculatedFvPatchScalarField::typeName
            )
        );

        PtrList<volScalarField> lawBiotCoeffs(laws.size());

        forAll(laws, lawI)
        {
            const poroMechanicalLaw2* poroLaw =
                dynamic_cast<const poroMechanicalLaw2*>(&laws[lawI]);
            const dimensionedScalar lawBiotCoeff
            (
                poroLaw
              ? poroLaw->biotCoeffValue()
              : dimensionedScalar("one", dimless, 1.0)
            );
            const tmp<volScalarField> tRho(laws[lawI].rho());
            const volScalarField& rho = tRho();

            lawBiotCoeffs.set
            (
                lawI,
                new volScalarField
                (
                    IOobject
                    (
                        "poroTractionLawBiotCoeff",
                        this->db().time().timeName(),
                        rho.db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    rho.mesh(),
                    lawBiotCoeff
                )
            );
        }

        mechanical.solSubMeshes().mapSubMeshVolFields<scalar>
        (
            lawBiotCoeffs,
            tb.ref()
        );

        return scalarField(tb().boundaryField()[patch().index()]);
    }

    // Update the coefficients associated with the patch field
    void poroTractionFvPatchVectorField::updateCoeffs()
    {

        if (updated()||this->db().time().timeIndex()<this->db().time().startTimeIndex()+1)
        {
            return;
        }

        traction() = tractionSeries_->value(this->db().time().timeOutputValue());
        pressure() = pressureSeries_->value(this->db().time().timeOutputValue());


        scalarField totalPressure(pressure());

        const solidModel &solMod = lookupSolidModel(patch().boundaryMesh().mesh());

        if(!totalTraction_)
        {
            const word pName(mechanicalPressureFieldName(solMod));
            const fvPatchField<scalar>& pPatch =
                this->patch().patchField<volScalarField, scalar>
                (
                    this->db().lookupObject<volScalarField>(pName)
                );
            const scalarField bPatch(patchBiotCoeff(solMod));

            totalPressure = pressure() + bPatch*pPatch;
        }
        pressure() = totalPressure;

        // Keep the inherited traction()/pressure() values in sync for solid
        // models that enforce solidTraction boundaries explicitly, but do not
        // call solidTraction::updateCoeffs(): it may refresh those values from
        // its own fields/series before computing the gradient.
        gradient() =
            relaxFac_ * solMod.tractionBoundarySnGrad(
                            traction(), pressure(), patch()) +
            (1.0 - relaxFac_) * gradient();

        fixedGradientFvPatchVectorField::updateCoeffs();
    }

    void poroTractionFvPatchVectorField::evaluate(
        const Pstream::commsTypes commsType)
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }

        // Lookup the gradient field
        const fvPatchField<tensor> &gradField =
            patch().lookupPatchField<volTensorField, tensor>(
                "grad(" + internalField().name() + ")"
            );

        // Face unit normals
        const vectorField n((patch().nf())());

        // Delta vectors
        const vectorField delta((patch().delta())());

        // Non-orthogonal correction vectors
        const vectorField k((((I - sqr(n)) & delta))());

        if (secondOrder_)
        {
            const vectorField dUP(((k & gradField.patchInternalField()))());
            const vectorField nGradUP(((n & gradField.patchInternalField()))());

            Field<vector>::operator=(
                patchInternalField() + dUP + 0.5 * (gradient() + nGradUP) / patch().deltaCoeffs());
        }
        else
        {

            Field<vector>::operator=(
                patchInternalField() + (k & gradField.patchInternalField()) + gradient() / patch().deltaCoeffs());
        }

        fvPatchField<vector>::evaluate();
    }

    void poroTractionFvPatchVectorField::write(Ostream &os) const
    {
        solidTractionFvPatchVectorField::write(os);

        os.writeKeyword("total")
            << totalTraction_ << token::END_STATEMENT << nl;
        os.writeKeyword("limitCoeff")
            << limitCoeff_ << token::END_STATEMENT << nl;
        if(!totalTraction_)
        {
            os.writeKeyword("pressureFieldName")
                << pressureFieldName_ << token::END_STATEMENT << nl;
        }

    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(fvPatchVectorField, poroTractionFvPatchVectorField);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
