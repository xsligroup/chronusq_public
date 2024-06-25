#pragma once

#include <chronusq_sys.hpp>
#include <fields.hpp>
#include <orbitalmodifieroptions.hpp>

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
  template < FieldEnvelopeType _Typ > struct FieldEnvelope;

  // FieldEnvelope typedefs
  using LinRampField   = FieldEnvelope<FieldEnvelopeType::LinRamp>;
  using GaussianField  = FieldEnvelope<FieldEnvelopeType::Gaussian>;
  using StepField      = FieldEnvelope<FieldEnvelopeType::Step>;
  using PlaneWaveField = FieldEnvelope<FieldEnvelopeType::PlaneWave>;

  /**
   *  \brief FieldEnvelope specification for a Linear Ramp function
   *  envelope.
   */ 
  template<>
  struct FieldEnvelope<FieldEnvelopeType::LinRamp> : FieldEnvelopeBase {

    FieldEnvelope()                      = delete;
    FieldEnvelope(const FieldEnvelope &) = default;
    FieldEnvelope(FieldEnvelope &&)      = default;

    FieldEnvelope(double on, double off): FieldEnvelopeBase(on,off){ };

    double getAmp(double t);

  }; // struct FieldEnvelop<LinRamp>


  /**
   *  \brief FieldEnvelop specification for a Gaussian function
   *  envelope.
   */ 
  template<>
  struct FieldEnvelope<FieldEnvelopeType::Gaussian> : FieldEnvelopeBase {

    FieldEnvelope()                      = delete;
    FieldEnvelope(const FieldEnvelope &) = default;
    FieldEnvelope(FieldEnvelope &&)      = default;

    FieldEnvelope(double on, double off): FieldEnvelopeBase(on,off){ };
    FieldEnvelope(double on, double off, double alpha): FieldEnvelopeBase(on,off), alpha(alpha) { };

    double getAmp(double t);

    double alpha; // Gaussian envelope parameter should be > 0. Since the negative sign is in getAmp.
    void setAlpha(double alpha) { this->alpha = alpha; };

  }; // struct FieldEnvelop<Gaussian>


  /**
   *  \brief FieldEnvelop specification for a step function
   *  envelope.
   */ 
  template<>
  struct FieldEnvelope<FieldEnvelopeType::Step> : FieldEnvelopeBase {

    FieldEnvelope()                      = delete;
    FieldEnvelope(const FieldEnvelope &) = default;
    FieldEnvelope(FieldEnvelope &&)      = default;

    FieldEnvelope(double on, double off): FieldEnvelopeBase(on,off){ };

    double getAmp(double t);

  }; // struct FieldEnvelop<Step>

  /**
   *  \brief FieldEnvelop specification for a plane wave function
   *  envelope.
   */ 
  template<>
  struct FieldEnvelope<FieldEnvelopeType::PlaneWave> : FieldEnvelopeBase {

    FieldEnvelope()                      = delete;
    FieldEnvelope(const FieldEnvelope &) = default;
    FieldEnvelope(FieldEnvelope &&)      = default;

    FieldEnvelope(double on, double off ): FieldEnvelopeBase(on,off) { };
    FieldEnvelope(double on, double off, double omega): FieldEnvelopeBase(on,off), omega(omega) { };
    FieldEnvelope(double on, double off, double omega, bool doCos): FieldEnvelopeBase(on,off), omega(omega), doCos(doCos) { };

    double getAmp(double t);

    double omega;
    bool doCos = true;
    void setOmega(double omega) { this->omega = omega; };
    void setdoCos(bool doCos) { this->doCos = doCos; };

  }; // struct FieldEnvelop<PlaneWave>


  
  /**
   *  \brief Cast a templated FieldEnvelope shared_ptr to
   *  Base class.
   *
   *  dynamic_pointer_cast FieldEnvelope<> -> FieldEnvelopBase
   *
   *  \param [in] x Shared pointer for a FieldEnvelope object
   *  \returns      Shared pointer for a FieldEnvelopeBase object
   */ 
  template < FieldEnvelopeType _Typ >
  std::shared_ptr<FieldEnvelopeBase> 
    cast(std::shared_ptr<FieldEnvelope<_Typ>> x) {
    return std::dynamic_pointer_cast<
             FieldEnvelopeBase,FieldEnvelope<_Typ>>(x);
  }; // cast(FieldEnvelope) -> FieldEnvelopeBase

}; // namespace ChronusQ

