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
    void yMinPoissonDirichlet_54215ac932abcdd2651c8025aff523b920a408cc(bool load)
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
    yMinPoissonDirichletFixedValueFvPatchScalarField
);


const char* const yMinPoissonDirichletFixedValueFvPatchScalarField::SHA1sum =
    "54215ac932abcdd2651c8025aff523b920a408cc";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

yMinPoissonDirichletFixedValueFvPatchScalarField::
yMinPoissonDirichletFixedValueFvPatchScalarField
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
        Info<<"construct yMinPoissonDirichlet sha1: 54215ac932abcdd2651c8025aff523b920a408cc"
            " from patch/dictionary\n";
    }
}


yMinPoissonDirichletFixedValueFvPatchScalarField::
yMinPoissonDirichletFixedValueFvPatchScalarField
(
    const yMinPoissonDirichletFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct yMinPoissonDirichlet sha1: 54215ac932abcdd2651c8025aff523b920a408cc"
            " from patch/DimensionedField/mapper\n";
    }
}


yMinPoissonDirichletFixedValueFvPatchScalarField::
yMinPoissonDirichletFixedValueFvPatchScalarField
(
    const yMinPoissonDirichletFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct yMinPoissonDirichlet sha1: 54215ac932abcdd2651c8025aff523b920a408cc "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

yMinPoissonDirichletFixedValueFvPatchScalarField::
~yMinPoissonDirichletFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy yMinPoissonDirichlet sha1: 54215ac932abcdd2651c8025aff523b920a408cc\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void yMinPoissonDirichletFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs yMinPoissonDirichlet sha1: 54215ac932abcdd2651c8025aff523b920a408cc\n";
    }

//{{{ begin code
    #line 177 "/home/a776/workdocuments/上交船舶/slover/student_project/cases/04_poisson_equation/01_poisson_manufactured_quad/N10/0/phi!boundaryField/yMin"
 
            const vectorField& Cf = patch().Cf();
            scalarField values(patch().size(), 0.0);
            forAll(Cf, i)
            {
                const scalar x = Cf[i].x();
                const scalar y = Cf[i].y();
                values[i] = cos(constant::mathematical::pi*x)
                          * cos(constant::mathematical::pi*y);
            }
            operator==(values);
        
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

