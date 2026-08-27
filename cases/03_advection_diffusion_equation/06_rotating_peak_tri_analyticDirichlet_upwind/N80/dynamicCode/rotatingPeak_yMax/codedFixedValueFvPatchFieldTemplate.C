/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "read.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // Unique function name that can be checked
    // to ensure the correct library version has been loaded
    void rotatingPeak_yMax_081206e55da32c49b81a6b735e0d1ecbf7f4525b(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    rotatingPeak_yMaxFixedValueFvPatchScalarField
);


const char* const rotatingPeak_yMaxFixedValueFvPatchScalarField::SHA1sum =
    "081206e55da32c49b81a6b735e0d1ecbf7f4525b";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatingPeak_yMaxFixedValueFvPatchScalarField::
rotatingPeak_yMaxFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct rotatingPeak_yMax sha1: 081206e55da32c49b81a6b735e0d1ecbf7f4525b"
            " from patch/dictionary\n";
    }
}


rotatingPeak_yMaxFixedValueFvPatchScalarField::
rotatingPeak_yMaxFixedValueFvPatchScalarField
(
    const rotatingPeak_yMaxFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct rotatingPeak_yMax sha1: 081206e55da32c49b81a6b735e0d1ecbf7f4525b"
            " from patch/DimensionedField/mapper\n";
    }
}


rotatingPeak_yMaxFixedValueFvPatchScalarField::
rotatingPeak_yMaxFixedValueFvPatchScalarField
(
    const rotatingPeak_yMaxFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct rotatingPeak_yMax sha1: 081206e55da32c49b81a6b735e0d1ecbf7f4525b "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rotatingPeak_yMaxFixedValueFvPatchScalarField::
~rotatingPeak_yMaxFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy rotatingPeak_yMax sha1: 081206e55da32c49b81a6b735e0d1ecbf7f4525b\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rotatingPeak_yMaxFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs rotatingPeak_yMax sha1: 081206e55da32c49b81a6b735e0d1ecbf7f4525b\n";
    }

//{{{ begin code
    #line 12989 "/home/a776/workdocuments/上交船舶/slover/student_project/cases/03_advection_diffusion_equation/06_rotating_peak_tri_analyticDirichlet_upwind/N80/0/phi!boundaryField/yMax"

            const scalar epsilon = 1.0000000000000000e-03;
            const scalar t0 = 1.5707963267948966e+00;
            const scalar x0 = 0.0000000000000000e+00;
            const scalar y0 = 5.0000000000000000e-01;
            const scalar tau = this->db().time().value();
            const scalar physicalTime = t0 + tau;
            const scalar xhat = x0*cos(tau) - y0*sin(tau);
            const scalar yhat = x0*sin(tau) + y0*cos(tau);
            const vectorField& faceCentres = this->patch().Cf();
            scalarField values(faceCentres.size(), Zero);
            forAll(faceCentres, faceI)
            {
                const scalar dx = faceCentres[faceI].x() - xhat;
                const scalar dy = faceCentres[faceI].y() - yhat;
                const scalar radius2 = dx*dx + dy*dy;
                values[faceI] =
                    exp(-radius2/(4.0*epsilon*physicalTime))
                   /(4.0*constant::mathematical::pi*epsilon*physicalTime);
            }
            operator==(values);
        
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

