/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#pragma once

#include <chronusq_sys.hpp>
#include <util/typedefs.hpp>

namespace ChronusQ {


  enum FieldGauge {
    Length,
    Velocity
  };

  enum EMFieldTyp {
    Electric,
    Magnetic
  };


  struct EMFieldBase {

    size_t     size;
    EMFieldTyp emFieldTyp;
    FieldGauge fieldGauge;

    EMFieldBase(EMFieldTyp em = Electric, FieldGauge fg = Length,
        size_t sz = 3): fieldGauge(fg), emFieldTyp(em), size(sz) { };

    virtual std::valarray<double> getAmp() = 0;

  };

  /**
   *  \brief Tensorial specification of an EM field perturbation.
   *
   *  Handles the specification of various descriptions of the EM field
   *  perturbation, such as the multipole expansion.
   *
   *  Expected field amplitudes of the following:
   *  1. Dipole operator     - [x, y, z]
   *  2. Quadrupole operator - [xx, xy, xz, yy, yz, zz]
   *  ect...
   */ 
  template <typename _VecTyp>
  struct EMField : virtual public EMFieldBase {

    _VecTyp ampVec;

    EMField()                 = delete;
    EMField( const EMField& ) = default;
    EMField( EMField&& )      = default;

    EMField(EMFieldTyp em, FieldGauge fg,
      const _VecTyp &x): ampVec(x), EMFieldBase{em,fg, x.size()}{ };

    EMField(EMFieldTyp em, const _VecTyp &x): 
      ampVec(x), EMFieldBase{em,Length, x.size()}{ };

    EMField(const _VecTyp &x): EMField(Electric,Length,x){ };

    std::valarray<double> getAmp() {
      return std::valarray<double>(ampVec.begin(),ampVec.size());
    }

  };


  // Useful aliases
  using DipoleField     = EMField<std::array<double,3>>;
  using QuadrupoleField = EMField<std::array<double,6>>;
  using OctupoleField   = EMField<std::array<double,10>>;








  /**
   *  \brief Cast a templated EMField shared_ptr to
   *  Base class.
   *
   *  dynamic_pointer_cast EMField<> -> EMFieldBase
   *
   *  \param [in] x Shared pointer for a EMField object
   *  \returns      Shared pointer for a EMFieldBase object
   */ 
  template < typename _UnitVector >
  std::shared_ptr<EMFieldBase> 
    cast(std::shared_ptr<EMField<_UnitVector>> x) {
    return std::dynamic_pointer_cast<
             EMFieldBase,EMField<_UnitVector>>(x);
  };
  
  /**
   *  \brief Full specification of a static EM perturbation
   *
   *  General to multiple fields.
   */ 
  struct EMPerturbation {

    std::vector<
      std::shared_ptr<EMFieldBase>
        > fields; ///< Fields for the EM perturbation


    inline void addField(std::shared_ptr<EMPerturbation> otherPert){
      for(auto &field : otherPert->fields) {
          this->addField(field);
      }
    }

    inline void addField(const EMPerturbation& otherPert) { addField(std::make_shared<EMPerturbation>(otherPert)); };

    /**
     *  \brief Add a field to the perturbation.
     *
     *  \param [in] f    Tensorial field specification
     */ 
    inline void addField(std::shared_ptr<EMFieldBase> f) {
      fields.push_back(f);
    };

    template <typename _Field>
    inline void addField(std::shared_ptr<_Field> f){ addField(cast(f)); };

    template <typename _Field>
    inline void addField(_Field f){ addField(std::make_shared<_Field>(f)); }

    template <typename _UV>
    inline void addField(EMFieldTyp em, FieldGauge fg, const _UV &uv) {
      addField(EMField<_UV>(em,fg,uv));    
    }

    template <typename _UV>
    inline void addField(EMFieldTyp em, const _UV &uv) {
      addField(EMField<_UV>(em,uv));    
    }



    template <size_t N>
    inline std::array<double, N> getNAmp(EMFieldTyp TYPE) {

      std::array<double, N> amp = {0., 0., 0.};


      for(auto &field : fields) {

        if( field->emFieldTyp == TYPE  and field->size == N ) {

          auto tmp_amp = field->getAmp();

          for( auto k = 0; k < N; k++) amp[k] += tmp_amp[k];

        }

      }

      return amp;

    }


    inline std::array<double, 3> getDipoleAmp(EMFieldTyp TYPE) {

      return getNAmp<3>( TYPE );

    };
    
  }; // struct EMPerturbation



  // Helper functions
  inline bool pert_has_type( const EMPerturbation &pert, 
      const EMFieldTyp typ ) {

    return std::any_of(pert.fields.begin(), pert.fields.end(),
        [&](std::shared_ptr<EMFieldBase> x) {
          return x->emFieldTyp == typ;
        });

  }


  // Instantiations of Amplitude getters
};

