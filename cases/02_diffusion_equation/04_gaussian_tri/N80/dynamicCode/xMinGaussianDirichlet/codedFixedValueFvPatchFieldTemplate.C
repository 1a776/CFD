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
    void xMinGaussianDirichlet_fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10(bool load)
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
    xMinGaussianDirichletFixedValueFvPatchScalarField
);


const char* const xMinGaussianDirichletFixedValueFvPatchScalarField::SHA1sum =
    "fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

xMinGaussianDirichletFixedValueFvPatchScalarField::
xMinGaussianDirichletFixedValueFvPatchScalarField
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
        Info<<"construct xMinGaussianDirichlet sha1: fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10"
            " from patch/dictionary\n";
    }
}


xMinGaussianDirichletFixedValueFvPatchScalarField::
xMinGaussianDirichletFixedValueFvPatchScalarField
(
    const xMinGaussianDirichletFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct xMinGaussianDirichlet sha1: fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10"
            " from patch/DimensionedField/mapper\n";
    }
}


xMinGaussianDirichletFixedValueFvPatchScalarField::
xMinGaussianDirichletFixedValueFvPatchScalarField
(
    const xMinGaussianDirichletFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct xMinGaussianDirichlet sha1: fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

xMinGaussianDirichletFixedValueFvPatchScalarField::
~xMinGaussianDirichletFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy xMinGaussianDirichlet sha1: fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void xMinGaussianDirichletFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs xMinGaussianDirichlet sha1: fa4a0ea4674e0a9b971fc0ecbe38b5362f204b10\n";
    }

//{{{ begin code
    #line 12829 "/home/a776/workdocuments/上交船舶/slover/student_project/cases/02_diffusion_equation/04_gaussian_tri/N80/0/phi!boundaryField/xMin"
 
            const scalar t = this->db().time().value();
            const scalar denom = 1.0 + 4.0*t;
            const scalar invDenom = 1.0/denom;
            const scalar mu = 1;
            const vectorField& Cf = patch().Cf();
            scalarField values(patch().size(), 0.0);
            forAll(Cf, i)
            {
                const scalar x = Cf[i].x();
                const scalar y = Cf[i].y();
                values[i] = invDenom*exp(-mu*(x*x + y*y)*invDenom);
            }
            operator==(values);
        
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

