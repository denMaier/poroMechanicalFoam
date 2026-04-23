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

Class
    Foam::varSatPoroHydraulicModel

Description
    The varSatPoroHydraulic model class extends the 
    poroHydraulic model class with saturation laws describing the connection
    between pressure, saturation and hydraulic conductivity in variably saturated
    porous media.

SourceFiles
    varSatPoroHydraulicModel.C

Author
    Denis Maier, BAW. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroHydraulicModel_H
#define varSatPoroHydraulicModel_H

#include "poroHydraulicModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
  /*---------------------------------------------------------------------------*\
                        Class varSatPoroHydraulicModel Declaration
\*---------------------------------------------------------------------------*/

  class varSatPoroHydraulicModel
      : public poroHydraulicModel
  {
    // Private data

       // const uniformDimensionedVectorField g_;

    //- Zone-wise saturation laws coupling pressure, saturation and kr.
    mutable PtrList<saturationLaw> saturationLawPtr_;

    //- Disallow default bitwise copy construct
    varSatPoroHydraulicModel(const varSatPoroHydraulicModel &);

    //- Disallow default bitwise assignment
    void operator=(const varSatPoroHydraulicModel &);

  public:
    //- Runtime type information
    TypeName("varSatPoroHydraulicModel");

    // Constructors
    //- Construct from an fvMesh
    varSatPoroHydraulicModel(
        const volScalarField &pField,
        const dimensionedVector& gravity);
    //varSatPoroHydraulicModel(
    //    const volScalarField &pField);

    //- Destructor
    ~varSatPoroHydraulicModel();

    // Member Functions

    // Access

    //- Return the zone-wise air-entry or reference pressure used by the SWCCs.
    virtual tmp<volScalarField> pStar() const;

    //- Return cell-effective conductivity k*kr(p).
    virtual const tmp<volScalarField> kEff(const volScalarField &p) const;

    //- Return face-effective conductivity kf*interpolate(kr(p)).
    virtual const tmp<surfaceScalarField> kEfff(const volScalarField &p) const;

    //- Return relative conductivity from the active saturation laws.
    virtual tmp<volScalarField> kr(const volScalarField &p) const;

    //- Return saturation from the active saturation laws.
    virtual tmp<volScalarField> S(const volScalarField &p) const;

    //- Return specific moisture capacity dS/dp from the active saturation laws.
    virtual tmp<volScalarField> C(const volScalarField &p) const;

    //- Write function must be defined for regIOobjects
    virtual bool writeData(Ostream&) const
    {
      return false;
    }

    virtual void writeCoeffFields() const
    {
      //saturationLawPtr_->write();
      //poroHydraulicModel::writeCoeffFields();
    }

  };

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
