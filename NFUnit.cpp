#include "NFUnit.h"
#include <vector>
#include <cmath>

// =============================================================================
// Five-unit NF Li/Co enrich-strip cascade, MOLE-FLOW decomposition.
//
// SCENARIO -> UNIT map (scenarioIndex+1 == unitIdx):
//   1 = NFF   (feed stage)
//   2 = NFE1  (enriching stage 1)
//   3 = NFE2  (enriching stage 2 -- Li product unit)
//   4 = NFS1  (stripping stage 1)
//   5 = NFS2  (stripping stage 2 -- Co product unit)
//
// STREAM legend (15 broken streams; each has molLi, molCo, molW):
//    1 = NFF inlet   (Mixer1: feed + s5 + s12)
//    2 = NFF ret     -> Mixer2 (s10, NFS1 inlet)
//    3 = NFF perm    -> Mixer3 (s4, NFE1 inlet)
//    4 = NFE1 inlet  (Mixer3: s3 + s8)
//    5 = NFE1 ret    -> recycle to Mixer1 (s1)
//    6 = NFE1 perm   -> NFE2 inlet (s7, direct)
//    7 = NFE2 inlet  (= s6)
//    8 = NFE2 ret    -> recycle to Mixer3 (s4)
//    9 = NFE2 perm   -> Li PRODUCT
//   10 = NFS1 inlet  (Mixer2: s2 + s15)
//   11 = NFS1 ret    -> NFS2 inlet (s13, direct)
//   12 = NFS1 perm   -> recycle to Mixer1 (s1)
//   13 = NFS2 inlet  (= s11)
//   14 = NFS2 ret    -> Co PRODUCT
//   15 = NFS2 perm   -> recycle to Mixer2 (s10)
//
//   Mixer1 (s1) :  s1  = feed + s5 + s12
//   Mixer3 (s4) :  s4  = s3 + s8
//   Mixer2 (s10):  s10 = s2 + s15
//   Pass   (s7) :  s7  = s6
//   Pass   (s13):  s13 = s11
//
// Unit -> (inlet stream, ret stream, perm stream):
//   1(NFF)  -> in=s1,  ret=s2,  perm=s3
//   2(NFE1) -> in=s4,  ret=s5,  perm=s6
//   3(NFE2) -> in=s7,  ret=s8,  perm=s9
//   4(NFS1) -> in=s10, ret=s11, perm=s12
//   5(NFS2) -> in=s13, ret=s14, perm=s15
// =============================================================================


NFUnit::NFUnit(BranchingStrategy branchingStrategy) : STModel() {
    this->branching_strategy = branchingStrategy;
    this->scenario_names = {
        ScenarioNames::SCENARIO1,
        ScenarioNames::SCENARIO2,
        ScenarioNames::SCENARIO3,
        ScenarioNames::SCENARIO4,
        ScenarioNames::SCENARIO5
    };
    this->scenario_name = ScenarioNames::SCENARIO1;
    this->probability = 1.0;

    // ------------------------------------------------------------------
    // FIRST-STAGE: 15 streams x (molLi, molCo, molW) = 45 variables.
    // Ordering per stream s (0-based block start = 3*(s-1)):
    //   +0 molLiStream{s}, +1 molCoStream{s}, +2 molWStream{s}
    // Bounds mirror flowsheet5NFMolFlow.py: Li/Co in (1e-6,1e4), W in (1e-3,5e6).
    // ------------------------------------------------------------------
    this->first_stage_IX.clear();
    for (int s = 1; s <= 15; ++s) {
        this->first_stage_IX.push_back(mc::Interval(1e-6, 1e4));   // molLiStream{s}
        this->first_stage_IX.push_back(mc::Interval(1e-6, 1e4));   // molCoStream{s}
        this->first_stage_IX.push_back(mc::Interval(1e-3, 5e6));   // molWStream{s}
    }

    // ------------------------------------------------------------------
    // SECOND-STAGE: unit-local NF variables (per scenario/unit).
    // Bounds mirror flowsheet5NFMolFlow.py second-stage Var declarations.
    // ------------------------------------------------------------------
    this->second_stage_IX = {
        mc::Interval(1e5, 40e5),     // dP              0  [Pa]
        mc::Interval(0.1, 200.0),    // memLength       1  [m]
        mc::Interval(0.95, 0.95),    // memWidth        2  [m] (pinned via equal bounds)
        mc::Interval(2e-4, 5e-3),    // channelHeight   3  [m]
        mc::Interval(1e-6, 100.0),   // flowOutRetND    4
        mc::Interval(1e-6, 100.0),   // flowOutPermND   5
        mc::Interval(1e-6, 100.0),   // cLiRetND        6
        mc::Interval(1e-6, 100.0),   // cCoRetND        7
        mc::Interval(1e-6, 100.0),   // cLiPermND       8
        mc::Interval(1e-6, 500.0),   // cCoPermND       9
        mc::Interval(1e-8, 1e-4),    // Jw             10  [m3/(m2 s)]
        mc::Interval(1e-4, 1.0),     // fLi            11
        mc::Interval(1e-4, 1.0),     // fCo            12
        mc::Interval(1e-6, 30000.0), // cLiIn          13  [mol/m3]
        mc::Interval(1e-6, 30000.0), // cCoIn          14  [mol/m3]
        mc::Interval(0.0, 1e9),      // capex          15
        mc::Interval(0.0, 1e9),      // opex1          16
        mc::Interval(0.0, 1e9)       // opex2          17
    };
}

void NFUnit::buildDAG() {
    // ---- shared constants (mirror flowsheet5NFMolFlow.py params) ----
    const double NUM_ELEMENTS   = 5000.0;
    const double SCALE_FAC_FLOW = 1.0 / 5000.0;
    const double SCALE_FAC_CONC = 50.0;

    const double rg              = 8.314;
    const double temp            = 273.15 + 25.0;
    const double densityH2O      = 997.0;                 // kg/m3
    const double densityH2OGm3   = 997.0 * 1000.0;        // g/m3
    const double viscosity       = 1e-3;                  // Pa.s
    const double molweightH2OgMol = 18.01528;             // g/mol
    const double molweightLi     = 6.941;                 // g/mol
    const double molweightCo     = 58.933;                // g/mol

    const double Lp              = 1.0e-11;               // m3/(m2 s Pa)
    const double sigmaLi         = 0.05;
    const double sigmaCo         = 0.85;
    const double pLi             = 5.0e-6;                // m/s
    const double pCo             = 2.0e-7;                // m/s
    const double nuLi            = 2.0;
    const double nuCo            = 3.0;
    const double diffLi          = 1.03e-9;               // m2/s
    const double diffCo          = 0.72e-9;               // m2/s

    const double costElec        = 0.08;                  // $/kWh
    const double costMembrane    = 50.0;                  // $/m2/yr
    const double pumpEff         = 0.80;
    const double daysOperation   = 350.0;
    const double cpMax           = 2.0;

    // ---- fresh feed (mole-flow basis) ----
    const double feedFlow   = 1.0;    // m3/s
    const double feedConcLi = 50.0;   // mol/m3
    const double feedConcCo = 100.0;  // mol/m3
    const double molLiFeed  = feedFlow * feedConcLi;
    const double molCoFeed  = feedFlow * feedConcCo;
    const double molWFeed   = feedFlow * (densityH2OGm3 - feedConcLi * molweightLi
                                          - feedConcCo * molweightCo) / molweightH2OgMol;

    for (size_t scenarioIndex = 0; scenarioIndex < this->scenario_names.size(); ++scenarioIndex) {
        auto scenarioName = this->scenario_names[scenarioIndex];
        const int firstStageCount  = static_cast<int>(this->first_stage_IX.size());
        const int secondStageCount = static_cast<int>(this->second_stage_IX.size());
        const int nvars = firstStageCount + secondStageCount;

        this->X[scenarioName].resize(nvars);
        for (int i = 0; i < nvars; ++i) {
            this->X[scenarioName][i].set(&this->DAG[scenarioName]);
        }
        auto& vars = this->X[scenarioName];

        // ---- first-stage handles: molLiStream[s], molCoStream[s], molWStream[s] ----
        // 1-based stream index -> block base = 3*(s-1)
        auto molLi = [&](int s) -> mc::FFVar& { return vars[3 * (s - 1) + 0]; };
        auto molCo = [&](int s) -> mc::FFVar& { return vars[3 * (s - 1) + 1]; };
        auto molW  = [&](int s) -> mc::FFVar& { return vars[3 * (s - 1) + 2]; };

        // ---- second-stage handles ----
        const int ss = firstStageCount;
        mc::FFVar& dP             = vars[ss + 0];
        mc::FFVar& memLength      = vars[ss + 1];
        mc::FFVar& memWidth       = vars[ss + 2];
        mc::FFVar& channelHeight  = vars[ss + 3];
        mc::FFVar& flowOutRetND   = vars[ss + 4];
        mc::FFVar& flowOutPermND  = vars[ss + 5];
        mc::FFVar& cLiRetND       = vars[ss + 6];
        mc::FFVar& cCoRetND       = vars[ss + 7];
        mc::FFVar& cLiPermND      = vars[ss + 8];
        mc::FFVar& cCoPermND      = vars[ss + 9];
        mc::FFVar& Jw             = vars[ss + 10];
        mc::FFVar& fLi            = vars[ss + 11];
        mc::FFVar& fCo            = vars[ss + 12];
        mc::FFVar& cLiIn          = vars[ss + 13];
        mc::FFVar& cCoIn          = vars[ss + 14];
        mc::FFVar& capex          = vars[ss + 15];
        mc::FFVar& opex1          = vars[ss + 16];
        mc::FFVar& opex2          = vars[ss + 17];

        const int unitIdx = static_cast<int>(scenarioIndex) + 1;

        // ---- pick this unit's inlet stream (mole flows) ----
        // 1->s1, 2->s4, 3->s7, 4->s10, 5->s13
        int inStream;
        if      (unitIdx == 1) inStream = 1;
        else if (unitIdx == 2) inStream = 4;
        else if (unitIdx == 3) inStream = 7;
        else if (unitIdx == 4) inStream = 10;
        else                    inStream = 13;

        mc::FFVar& molLiIn = molLi(inStream);
        mc::FFVar& molCoIn = molCo(inStream);
        mc::FFVar& molWIn  = molW(inStream);

        // mol/s -> (flow [m3/s], conc [mol/m3]); flowIn divides by CONSTANT density
        mc::FFVar massIn = molLiIn * molweightLi + molCoIn * molweightCo + molWIn * molweightH2OgMol; // g/s
        mc::FFVar flowIn = massIn / densityH2OGm3;                                                    // m3/s
        mc::FFVar flowInPerElem = flowIn / NUM_ELEMENTS;

        // ---- constraint containers (same pattern as ED example) ----
        std::vector<mc::FFVar> constraints;
        auto addEqualityConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
            constraints.push_back(-expr);
        };
        auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
        };

        // ---- geometry ----
        mc::FFVar memArea      = memLength * memWidth;
        mc::FFVar memAreaTotal = memArea * NUM_ELEMENTS;

        // ---- hydrodynamics ----
        mc::FFVar avgFlowRet = 0.5 * (flowOutRetND * SCALE_FAC_FLOW + flowInPerElem);
        mc::FFVar uFeed      = avgFlowRet / (memWidth * channelHeight);
        mc::FFVar Re         = densityH2O * uFeed * (2.0 * channelHeight) / viscosity;

        const double ScLi = viscosity / (densityH2O * diffLi);
        const double ScCo = viscosity / (densityH2O * diffCo);

        mc::FFVar ShLi = 1.85 * pow(Re * ScLi * 2.0 * channelHeight / memLength, 1.0 / 3.0);
        mc::FFVar ShCo = 1.85 * pow(Re * ScCo * 2.0 * channelHeight / memLength, 1.0 / 3.0);

        mc::FFVar kLi = ShLi * diffLi / (2.0 * channelHeight);
        mc::FFVar kCo = ShCo * diffCo / (2.0 * channelHeight);

        // ---- concentration polarization ----
        mc::FFVar cpFactorLi = exp(Jw / kLi);
        mc::FFVar cpFactorCo = exp(Jw / kCo);

        mc::FFVar avgConcLiRet = 0.5 * (cLiRetND * SCALE_FAC_CONC + cLiIn);
        mc::FFVar avgConcCoRet = 0.5 * (cCoRetND * SCALE_FAC_CONC + cCoIn);

        mc::FFVar cLiMem = avgConcLiRet * cpFactorLi;
        mc::FFVar cCoMem = avgConcCoRet * cpFactorCo;

        // ---- osmotic pressure (sigma-weighted, effective) ----
        mc::FFVar piRetEff  = rg * temp * (sigmaLi * nuLi * avgConcLiRet + sigmaCo * nuCo * avgConcCoRet);
        mc::FFVar piPermEff = rg * temp * (sigmaLi * nuLi * cLiPermND * SCALE_FAC_CONC
                                           + sigmaCo * nuCo * cCoPermND * SCALE_FAC_CONC);

        // ---- limiting flux ----
        mc::FFVar JwCritCo = kCo * log(cpMax);

        // ---- feed-channel pressure drop ----
        mc::FFVar fFeed = 96.0 / Re;
        mc::FFVar presDropChannel = densityH2O * fFeed * memLength * pow(uFeed, 2.0) / (4.0 * channelHeight);

        // ---- cost intermediates ----
        mc::FFVar capexMembrane = costMembrane * memAreaTotal;
        mc::FFVar capexHousing  = 20.0 * memAreaTotal;
        mc::FFVar pumpPowerTMP     = dP * flowIn / pumpEff;
        mc::FFVar pumpPowerChannel = presDropChannel * avgFlowRet * NUM_ELEMENTS / pumpEff;

        // =================================================================
        // FIRST-STAGE mixers / passthroughs / product targets.
        // Added in EVERY scenario (not gated by unitIdx) -- matching the ED
        // example, where the flowsheet coupling constraints are contributed to
        // each scenario's subproblem DAG so every subproblem sees the full
        // first-stage feasible region. All operands here are first-stage vars.
        // =================================================================
        // Mixer1 (s1): s1 = feed + s5 + s12
        addEqualityConstraint(molLi(1) - molLiFeed - molLi(5) - molLi(12));
        addEqualityConstraint(molCo(1) - molCoFeed - molCo(5) - molCo(12));
        addEqualityConstraint(molW(1)  - molWFeed  - molW(5)  - molW(12));

        // Mixer3 (s4): s4 = s3 + s8
        addEqualityConstraint(molLi(4) - molLi(3) - molLi(8));
        addEqualityConstraint(molCo(4) - molCo(3) - molCo(8));
        addEqualityConstraint(molW(4)  - molW(3)  - molW(8));

        // Mixer2 (s10): s10 = s2 + s15
        addEqualityConstraint(molLi(10) - molLi(2) - molLi(15));
        addEqualityConstraint(molCo(10) - molCo(2) - molCo(15));
        addEqualityConstraint(molW(10)  - molW(2)  - molW(15));

        // Pass (s7 = s6)
        addEqualityConstraint(molLi(7) - molLi(6));
        addEqualityConstraint(molCo(7) - molCo(6));
        addEqualityConstraint(molW(7)  - molW(6));

        // Pass (s13 = s11)
        addEqualityConstraint(molLi(13) - molLi(11));
        addEqualityConstraint(molCo(13) - molCo(11));
        addEqualityConstraint(molW(13)  - molW(11));

        // Product purity / recovery targets.
        //   Li product = s9 (NFE2 perm); Co product = s14 (NFS2 ret).
        //   purity molLi9 >= 0.95*(molLi9+molCo9)  ->  <= form:
        //          0.95*(molLi9+molCo9) - molLi9 <= 0
        addLessEqualConstraint(0.95 * (molLi(9) + molCo(9)) - molLi(9));       // Li purity >= 95%
        addLessEqualConstraint(0.95 * (molLi(14) + molCo(14)) - molCo(14));    // Co purity >= 95%
        addLessEqualConstraint(0.95 * molLiFeed - molLi(9));                    // Li recovery >= 95%
        addLessEqualConstraint(0.95 * molCoFeed - molCo(14));                   // Co recovery >= 95%

        // =================================================================
        // UNIT PHYSICS (constraint1..constraint10 + bilinear cLiIn/cCoIn + cost)
        // ported unchanged from flowsheet5NFMolFlow.py
        // =================================================================

        // bilinear inlet-concentration definitions: cIn * flowIn == molIn
        addEqualityConstraint(cLiIn * flowIn - molLiIn);
        addEqualityConstraint(cCoIn * flowIn - molCoIn);

        // constraint1: volume balance (per element, nondimensional)
        addEqualityConstraint(flowOutPermND + flowOutRetND - flowInPerElem / SCALE_FAC_FLOW);

        // constraint2: Li component balance
        addEqualityConstraint(flowOutPermND * cLiPermND + flowOutRetND * cLiRetND
                              - (flowInPerElem * cLiIn) / (SCALE_FAC_FLOW * SCALE_FAC_CONC));

        // constraint3: Co component balance
        addEqualityConstraint(flowOutPermND * cCoPermND + flowOutRetND * cCoRetND
                              - (flowInPerElem * cCoIn) / (SCALE_FAC_FLOW * SCALE_FAC_CONC));

        // constraint4: solvent flux (Darcy w/ effective osmotic pressure)
        addEqualityConstraint(Jw - Lp * (dP - (piRetEff - piPermEff)));

        // constraint5: permeate volume = Jw * area
        addEqualityConstraint(flowOutPermND * SCALE_FAC_FLOW - Jw * memArea);

        // constraint6/7: Spiegler-Kedem F factors
        addEqualityConstraint(fLi - exp(-Jw * (1.0 - sigmaLi) / pLi));
        addEqualityConstraint(fCo - exp(-Jw * (1.0 - sigmaCo) / pCo));

        // constraint8/9: permeate concentration (Spiegler-Kedem transport)
        addEqualityConstraint(cLiPermND * SCALE_FAC_CONC * (1.0 - sigmaLi * fLi)
                              - (1.0 - sigmaLi) * cLiMem);
        addEqualityConstraint(cCoPermND * SCALE_FAC_CONC * (1.0 - sigmaCo * fCo)
                              - (1.0 - sigmaCo) * cCoMem);

        // constraint10: limiting-flux inequality  Jw/JwCritCo - 1 <= 0
        addLessEqualConstraint(Jw / JwCritCo - 1.0);

        // cost defining equalities
        addEqualityConstraint(capex - (capexMembrane + capexHousing));
        addEqualityConstraint(opex1 - costElec * pumpPowerTMP * (24.0 * daysOperation) / 1000.0);
        addEqualityConstraint(opex2 - costElec * pumpPowerChannel * (24.0 * daysOperation) / 1000.0);

        // =================================================================
        // LINKING: broken-stream mole flows == this unit's outlet mole flows.
        // Ret/Perm stream mole flows from ND outlet vars:
        //   flowRet  = flowOutRetND  * SCALE_FAC_FLOW * NUM_ELEMENTS
        //   cLiRet   = cLiRetND      * SCALE_FAC_CONC   (etc.)
        //   molLiRet = cLiRet  * flowRet
        //   molWRet  = flowRet * (densityH2OGm3 - cLiRet*molweightLi - cCoRet*molweightCo)/molweightH2OgMol
        // =================================================================
        mc::FFVar flowRet  = flowOutRetND  * SCALE_FAC_FLOW * NUM_ELEMENTS;
        mc::FFVar flowPerm = flowOutPermND * SCALE_FAC_FLOW * NUM_ELEMENTS;
        mc::FFVar cLiRet   = cLiRetND  * SCALE_FAC_CONC;
        mc::FFVar cCoRet   = cCoRetND  * SCALE_FAC_CONC;
        mc::FFVar cLiPerm  = cLiPermND * SCALE_FAC_CONC;
        mc::FFVar cCoPerm  = cCoPermND * SCALE_FAC_CONC;

        mc::FFVar molLiRet  = cLiRet * flowRet;
        mc::FFVar molCoRet  = cCoRet * flowRet;
        mc::FFVar molWRet   = flowRet * (densityH2OGm3 - cLiRet * molweightLi - cCoRet * molweightCo) / molweightH2OgMol;
        mc::FFVar molLiPerm = cLiPerm * flowPerm;
        mc::FFVar molCoPerm = cCoPerm * flowPerm;
        mc::FFVar molWPerm  = flowPerm * (densityH2OGm3 - cLiPerm * molweightLi - cCoPerm * molweightCo) / molweightH2OgMol;

        // unit -> (ret stream, perm stream)
        int retStream, permStream;
        if      (unitIdx == 1) { retStream = 2;  permStream = 3;  }
        else if (unitIdx == 2) { retStream = 5;  permStream = 6;  }
        else if (unitIdx == 3) { retStream = 8;  permStream = 9;  }
        else if (unitIdx == 4) { retStream = 11; permStream = 12; }
        else                    { retStream = 14; permStream = 15; }

        addEqualityConstraint(molLi(retStream)  - molLiRet);
        addEqualityConstraint(molCo(retStream)  - molCoRet);
        addEqualityConstraint(molW(retStream)   - molWRet);
        addEqualityConstraint(molLi(permStream) - molLiPerm);
        addEqualityConstraint(molCo(permStream) - molCoPerm);
        addEqualityConstraint(molW(permStream)  - molWPerm);

        // ---- objective: per-scenario unit cost (matches ED buildDAG: fresh per
        //      scenario, NOT accumulated across scenarios) ----
        mc::FFVar objective = this->probability * (capex + opex1 + opex2);

        // ---- attach this scenario's functions (objective first, then cons) ----
        std::vector<mc::FFVar> functions;
        functions.push_back(objective);
        functions.insert(functions.end(), constraints.begin(), constraints.end());
        this->F[scenarioName] = functions;
    }

    this->full_model_built = true;
}
void NFUnit::buildFullModelDAG(){

};
Ipopt::SmartPtr<STModel> NFUnit::clone() {
    Ipopt::SmartPtr<NFUnit> copy = new NFUnit();

    copy->scenario_name    = this->scenario_name;
    copy->first_stage_IX   = this->first_stage_IX;
    copy->second_stage_IX  = this->second_stage_IX;
    copy->scenario_names   = this->scenario_names;
    copy->probability      = this->probability;
    copy->clearDAG();
    if (this->full_model_built) {
        copy->buildFullModelDAG();
    } else {
        copy->buildDAG();
    }
    return copy;
}