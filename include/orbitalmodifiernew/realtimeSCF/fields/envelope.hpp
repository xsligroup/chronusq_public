#pragma once

#include <chronusq_sys.hpp>
#include <fields.hpp>
#include <realtime/enums.hpp>

namespace ChronusQ {

  /**
   *  \brief Base struct for the specification of a field envelope.
   *
   *  Provides a minimal structure to build a FieldEnvelope specification.
   *  All FieldEnvelope classes are derived from this one.
   */ 
  struct FieldEnvelopeBase {

    double tOn;  ///< Time to turn on the perturbation
    double tOff; ///< Time to turn off the perturbation

    // Default / Delete constructors
    FieldEnvelopeBase()                          = delete;
    FieldEnvelopeBase(const FieldEnvelopeBase &) = default;
    FieldEnvelopeBase(FieldEnvelopeBase &&)      = default;

    /**
     *  \brief FieldEnvelopeBase constructor.
     *
     *  \param [in] on   Populates tOn
     *  \param [in] off  Populated tOff
     */ 
    FieldEnvelopeBase(double on, double off): tOn(on), tOff(off){ };




    /**
     *  \brief Obtain the field amplitude at a specified time.
     *
     *  Pure virtual, to be specified by derived classes
     *
     *  \param [in] t Time to evaluate the field amplitude
     *  \returns      The scalar field amplitude at time \f$ t \f$.
     */ 
    virtual double getAmp(double t) = 0; 


  }; // struct FieldEnvelopeBase


  // Forward declaration of FieldEnvelope templates
  template < FieldEnvelopeTyp _Typ > struct FieldEnvelope;

  // FieldEnvelope typedefs
  using ConstantField  = FieldEnvelope<Constant>;
  using LinRampField   = FieldEnvelope<LinRamp>;
  using GaussianField  = FieldEnvelope<Gaussian>;
  using StepField      = FieldEnvelope<Step>;



  /**
   *  \brief FieldEnvelop specification for a step function
   *  envelope.
   */ 
  template<>
  struct FieldEnvelope<Step> : FieldEnvelopeBase {

    FieldEnvelope()                      = delete;
    FieldEnvelope(const FieldEnvelope &) = default;
    FieldEnvelope(FieldEnvelope &&)      = default;

    FieldEnvelope(double on, double off): FieldEnvelopeBase(on,off){ };

    double getAmp(double t);

  }; // struct FieldEnvelop<Step>

  
  /**
   *  \brief Cast a templated FieldEnvelope shared_ptr to
   *  Base class.
   *
   *  dynamic_pointer_cast FieldEnvelope<> -> FieldEnvelopBase
   *
   *  \param [in] x Shared pointer for a FieldEnvelope object
   *  \returns      Shared pointer for a FieldEnvelopeBase object
   */ 
  template < FieldEnvelopeTyp _Typ >
  std::shared_ptr<FieldEnvelopeBase> 
    cast(std::shared_ptr<FieldEnvelope<_Typ>> x) {
    return std::dynamic_pointer_cast<
             FieldEnvelopeBase,FieldEnvelope<_Typ>>(x);
  }; // cast(FieldEnvelope) -> FieldEnvelopeBase

}; // namespace ChronusQ

