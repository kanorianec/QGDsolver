/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver library, based on OpenFOAM+.

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

Class
    Foam::qgd::RSWEModel

Group
    grpQGDCoeffs

Description
    Class for one of possible ways of $\tau$ calculating.
    in this case $\tau$ parameter is calculating by

    $$ \tau = \alpha \cdot \frac{h}{U_0} $$

    Where:

    $\tau$ is regularisation parameter, $\alpha$ is parameter
    which depends by user in constant foulder, $h$ is size of the mesh
    and $U_0$ is characteristic velocity

SourceFiles
    RSWEModel.C

\*---------------------------------------------------------------------------*/

#include "QGDCoeffs.H"

#ifndef RSWEModel_H
#define RSWEModel_H

namespace Foam
{

namespace qgd
{

class RSWEModel : public QGDCoeffs
{

public:

    //-
    TypeName ("RSWEModel");

    //-
    RSWEModel
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~RSWEModel();

    //-
    void correct(const QGDThermo&);

private:

    //-

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
