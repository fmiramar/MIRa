// SPDX-License-Identifier: GPL-3.0-or-later
#include "MIR2Plugin.hpp"
#include "mir2/Harmonic.hpp"
namespace {
using namespace mir2plugin;
inline int boundedInt(double value,int low,int high){return mir2::finite(value)?int(mir2::clamp(value,double(low),double(high))):low;}
enum class HarmonicKind{peaks,inharmonicity,tristimulus,oddEven,energy,centroid,deviation,spread,variation,noisiness,roughness};
struct HarmonicUnit:Unit{double* magnitude;int fftSize,capacity;mir2::SpectralPeak peaks[128];mir2::HarmonicFrame previous;bool ready;};
void harmonicDestroy(HarmonicUnit* unit){release(unit->mWorld,unit->magnitude);}
template<HarmonicKind kind>void harmonicNext(HarmonicUnit* unit,int){
    if(!(IN0(0)>=0))return;
    SndBuf* buffer=resolveBuffer(unit,0);if(!buffer)return;
    {LOCK_SNDBUF(buffer);if(!fillMagnitudes(buffer,unit->fftSize,unit->magnitude))return;}
    const double binHz=unit->mWorld->mSampleRate/unit->fftSize;
    const int count=mir2::extractSpectralPeaks(unit->magnitude,unit->fftSize/2+1,binHz,20,0,
        mir2::sanitize(IN0(4)),unit->peaks,kind==HarmonicKind::roughness?unit->capacity:128);
    if constexpr (kind==HarmonicKind::roughness) {
        ZOUT0(0)=float(mir2::sensoryRoughness(unit->peaks,count,boundedInt(IN0(1),0,1),IN0(3)>0));
    } else {
        const auto frame=mir2::associateHarmonics(unit->peaks,count,mir2::sanitize(IN0(1)),unit->capacity,IN0(3));
        if constexpr (kind==HarmonicKind::peaks) {
            for(int i=0;i<unit->capacity;++i){ZOUT0(i)=float(frame.frequencies[i]);ZOUT0(unit->capacity+i)=float(frame.amplitudes[i]);}
        } else {
            const auto summary=mir2::summarizeHarmonics(frame,unit->ready?&unit->previous:nullptr);
            if constexpr (kind==HarmonicKind::tristimulus) {
                for(int i=0;i<3;++i)ZOUT0(i)=float(summary.tristimulus[i]);
            } else {
                double result=0;
                if constexpr (kind==HarmonicKind::inharmonicity) result=summary.inharmonicity;
                else if constexpr (kind==HarmonicKind::oddEven) result=summary.oddEven;
                else if constexpr (kind==HarmonicKind::centroid) result=summary.centroid;
                else if constexpr (kind==HarmonicKind::deviation) result=summary.deviation;
                else if constexpr (kind==HarmonicKind::spread) result=summary.spread;
                else if constexpr (kind==HarmonicKind::variation) {
                    result=summary.variation;unit->previous=frame;unit->ready=frame.count>0;
                } else if constexpr (kind==HarmonicKind::energy||kind==HarmonicKind::noisiness) {
                    result=mir2::harmonicEnergyRatio(unit->magnitude,unit->fftSize/2+1,frame);
                    if constexpr (kind==HarmonicKind::noisiness) {
                        double energy=0;
                        for(int i=0;i<=unit->fftSize/2;++i)energy+=unit->magnitude[i]*unit->magnitude[i];
                        result=energy>1e-20?1-result:0;
                    }
                }
                ZOUT0(0)=float(mir2::sanitize(result));
            }
        }
    }
}

template<HarmonicKind kind>void harmonicInit(HarmonicUnit* unit){
    unit->magnitude=nullptr;unit->fftSize=boundedInt(IN0(5),0,8192);unit->capacity=boundedInt(IN0(2),1,kind==HarmonicKind::roughness?128:64);
    unit->previous={};unit->ready=false;
    const int outputs=kind==HarmonicKind::peaks?2*unit->capacity:kind==HarmonicKind::tristimulus?3:1;
    if(mir2::isPowerOfTwo(unit->fftSize)&&unit->fftSize>=256&&int(unit->mNumOutputs)==outputs)
        unit->magnitude=allocate<double>(unit->mWorld,static_cast<std::size_t>(unit->fftSize/2+1));
    if(!unit->magnitude){Print("MIR2 harmonic: invalid dimensions or insufficient RT memory\n");SETCALC(ClearUnitOutputs);unit->mDone=true;}
    else SETCALC(harmonicNext<kind>);
    ClearUnitOutputs(unit,1);
}
#define HARMONIC_UNIT(Name,Kind)struct Name:HarmonicUnit{};void Name##_Ctor(Name* u){harmonicInit<HarmonicKind::Kind>(u);}void Name##_Dtor(Name* u){harmonicDestroy(u);}
HARMONIC_UNIT(HarmonicPeaks,peaks) HARMONIC_UNIT(Inharmonicity,inharmonicity) HARMONIC_UNIT(Tristimulus,tristimulus)
HARMONIC_UNIT(OddEvenRatio,oddEven) HARMONIC_UNIT(HarmonicEnergyRatio,energy) HARMONIC_UNIT(HarmonicCentroid,centroid)
HARMONIC_UNIT(HarmonicDeviation,deviation) HARMONIC_UNIT(HarmonicSpread,spread) HARMONIC_UNIT(HarmonicVariation,variation)
HARMONIC_UNIT(Noisiness,noisiness) HARMONIC_UNIT(SensoryDissonance2,roughness)
#undef HARMONIC_UNIT
}
namespace mir2plugin{void registerHarmonicUnits(){
    DefineDtorUnit(HarmonicPeaks);DefineDtorUnit(Inharmonicity);DefineDtorUnit(Tristimulus);DefineDtorUnit(OddEvenRatio);
    DefineDtorUnit(HarmonicEnergyRatio);DefineDtorUnit(HarmonicCentroid);DefineDtorUnit(HarmonicDeviation);
    DefineDtorUnit(HarmonicSpread);DefineDtorUnit(HarmonicVariation);DefineDtorUnit(Noisiness);DefineDtorUnit(SensoryDissonance2);
}}
