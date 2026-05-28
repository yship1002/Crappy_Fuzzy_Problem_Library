#include "EDUnits.h"

#include <vector>
#include <cmath>

EDUnits::EDUnits(BranchingStrategy branchingStrategy) : STModel() {
    this->branching_strategy = branchingStrategy;
    this->scenario_names = {
        ScenarioNames::SCENARIO1,
        ScenarioNames::SCENARIO2,
        ScenarioNames::SCENARIO3
    };
    this->scenario_name = ScenarioNames::SCENARIO1;
    this->probability = 1.0;

    this->first_stage_IX = {
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500),
        mc::Interval(1e-2, 10),
        mc::Interval(1e-1, 500)
    };

    this->second_stage_IX = {
        mc::Interval(1e-1, 1e3),
        mc::Interval(0, 1),
        mc::Interval(0.01, 100),
        mc::Interval(0.01, 10),
        mc::Interval(0.001, 1),
        mc::Interval(0.001, 1),
        mc::Interval(1e-4, 100),
        mc::Interval(1e-4, 100),
        mc::Interval(1e-4, 100),
        mc::Interval(1e-4, 100),
        mc::Interval(0, 5),
        mc::Interval(1e-8, 1),
        mc::Interval(1e-8, 1),
        mc::Interval(1e-7, 1),
        mc::Interval(1e-7, 1),
        mc::Interval(0, 1e6),
        mc::Interval(0, 1e6),
        mc::Interval(0, 1e6),
        mc::Interval(0, 1e6)
    };
}

void EDUnits::buildDAG() {
    const double numCells = 500.0;
    const double daysOperation = 7.0 * 365.0;
    const double rg = 8.314;
    const double temp = 273.15 + 25.0;
    const double faraday = 96485.3321;
    const double molweightN = 14.0067;
    const double molweightNh4Cl = 53.491;
    const double thicknessCem = 0.00011;
    const double thicknessAem = 0.000105;
    const double permSelCem = 0.9;
    const double permSelAem = 0.9;
    const double transCem = 0.95;
    const double transAem = 0.95;
    const double transIonConc = 0.491;
    const double transIonDil = 0.509;
    const double saltDiffCem = 1e-12;
    const double saltDiffAem = 1e-12;
    const double saltDiffConc = 1.84e-9;
    const double saltDiffDil = 1.84e-9;
    const double waterPermCem = 7.79e-11 / 3600.0;
    const double waterPermAem = 6.29e-11 / 3600.0;
    const double vantHoffNumber = 2.0;
    const double molweightH2O = 0.018;
    const double molweightH20gMol = 18.01528;
    const double densityH2O = 1000.0;
    const double densityH2OGm3 = 1000.0 * 1000.0;
    const double densityManure = 1000.0;
    const double viscosityManure = 1e-3;
    const double waterTransNumber = 6.0 + 8.0 + 2.0;
    const double osmoticCoeff = 1.0;
    const double resBlank = 2.5e-5;
    const double resCem = 2.5e-4;
    const double resAem = 1.8e-4;
    const double conducConcentrate = 149.6e-4;
    const double conducDilute = 149.6e-4;
    const double costElec = 0.08147;

    const double flowFeed = 192555.0 / (1000.0 * 24.0 * 60.0 * 60.0);
    const double concFeed = 1366.0;
    const double molNFeed = flowFeed * concFeed / molweightN;
    const double molWFeed = flowFeed * (densityH2OGm3 - concFeed) / molweightH20gMol;

    const double scaleFacFlowByUnit[4] = {
        0.0,
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0),
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0),
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0)
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        1366.0 / molweightN,
        1366.0 / molweightN,
        1366.0 / molweightN
    };

    for (size_t scenarioIndex = 0; scenarioIndex < this->scenario_names.size(); ++scenarioIndex) {
        auto scenarioName = this->scenario_names[scenarioIndex];
        const int firstStageCount = static_cast<int>(this->first_stage_IX.size());
        const int secondStageCount = static_cast<int>(this->second_stage_IX.size());
        const int nvars = firstStageCount + secondStageCount;

        this->X[scenarioName].resize(nvars);
        for (int i = 0; i < nvars; ++i) {
            this->X[scenarioName][i].set(&this->DAG[scenarioName]);
        }

        auto& vars = this->X[scenarioName];

        // Variable ordering: molNStream1, molWStream1, ..., molNStream8, molWStream8
        mc::FFVar& molNStream1 = vars[0];
        mc::FFVar& molWStream1 = vars[1];
        mc::FFVar& molNStream2 = vars[2];
        mc::FFVar& molWStream2 = vars[3];
        mc::FFVar& molNStream3 = vars[4];
        mc::FFVar& molWStream3 = vars[5];
        mc::FFVar& molNStream4 = vars[6];
        mc::FFVar& molWStream4 = vars[7];
        mc::FFVar& molNStream5 = vars[8];
        mc::FFVar& molWStream5 = vars[9];
        mc::FFVar& molNStream6 = vars[10];
        mc::FFVar& molWStream6 = vars[11];
        mc::FFVar& molNStream7 = vars[12];
        mc::FFVar& molWStream7 = vars[13];
        mc::FFVar& molNStream8 = vars[14];
        mc::FFVar& molWStream8 = vars[15];

        const int secondStageStart = firstStageCount;

        mc::FFVar& current = vars[secondStageStart + 0];
        mc::FFVar& flowSplit = vars[secondStageStart + 1];
        mc::FFVar& memLength = vars[secondStageStart + 2];
        mc::FFVar& memWidth = vars[secondStageStart + 3];
        mc::FFVar& thicknessConcentrate = vars[secondStageStart + 4];
        mc::FFVar& thicknessDilute = vars[secondStageStart + 5];
        mc::FFVar& flowOutConcentrateNd = vars[secondStageStart + 6];
        mc::FFVar& flowOutDiluteNd = vars[secondStageStart + 7];
        mc::FFVar& concOutConcentrateNd = vars[secondStageStart + 8];
        mc::FFVar& concOutDiluteNd = vars[secondStageStart + 9];
        mc::FFVar& voltCellPair = vars[secondStageStart + 10];
        mc::FFVar& resConcentrate = vars[secondStageStart + 11];
        mc::FFVar& resDilute = vars[secondStageStart + 12];
        mc::FFVar& voltNonOhmicCem = vars[secondStageStart + 13];
        mc::FFVar& voltNonOhmicAem = vars[secondStageStart + 14];
        mc::FFVar& capex = vars[secondStageStart + 15];
        mc::FFVar& opex1 = vars[secondStageStart + 16];
        mc::FFVar& opex2 = vars[secondStageStart + 17];
        mc::FFVar& costTotal = vars[secondStageStart + 18];

        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar massStream1 = molNStream1 * molweightN + molWStream1 * molweightH20gMol;
        mc::FFVar flowStream1 = massStream1 / densityH2OGm3;
        mc::FFVar concStream1 = molNStream1 / flowStream1;

        mc::FFVar massStream4 = molNStream4 * molweightN + molWStream4 * molweightH20gMol;
        mc::FFVar flowStream4 = massStream4 / densityH2OGm3;
        mc::FFVar concStream4 = molNStream4 / flowStream4;

        mc::FFVar massStream5 = molNStream5 * molweightN + molWStream5 * molweightH20gMol;
        mc::FFVar flowStream5In = massStream5 / densityH2OGm3;
        mc::FFVar concStream5In = molNStream5 / flowStream5In;

        mc::FFVar flowInEd = flowStream1;
        mc::FFVar concInEd = concStream1;
        if (unitIdx == 2) {
            flowInEd = flowStream4;
            concInEd = concStream4;
        } else if (unitIdx == 3) {
            flowInEd = flowStream5In;
            concInEd = concStream5In;
        }

        mc::FFVar flowInConc = flowSplit * flowInEd / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInEd / numCells;
        mc::FFVar concInConc = concInEd;
        mc::FFVar concInDil = concInEd;

        mc::FFVar uConc = 0.5 * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * thicknessConcentrate);
        mc::FFVar uDil = 0.5 * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * thicknessDilute);

        mc::FFVar reConc = densityManure * uConc * (2.0 * thicknessConcentrate) / viscosityManure;
        mc::FFVar reDil = densityManure * uDil * (2.0 * thicknessDilute) / viscosityManure;

        mc::FFVar scConc = viscosityManure / (densityManure * saltDiffConc);
        mc::FFVar scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConcCem = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDilCem = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);
        mc::FFVar shConcAem = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDilAem = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInConc);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInDil);

        mc::FFVar iLimConcCem = shConcCem * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transCem - transIonConc));
        mc::FFVar iLimDilCem = shDilCem * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transCem - transIonConc));
        mc::FFVar iLimConcAem = shConcAem * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transAem - transIonDil));
        mc::FFVar iLimDilAem = shDilAem * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transAem - transIonDil));

        mc::FFVar concInConcIntCem = concInConc * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concInDilIntCem = concInDil * (1.0 - current / (memLength * memWidth) / iLimDilCem);
        mc::FFVar concInConcIntAem = concInConc * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        mc::FFVar concInDilIntAem = concInDil * (1.0 - current / (memLength * memWidth) / iLimDilAem);

        mc::FFVar concOutConcIntCem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concOutDilIntCem = concOutDiluteNd * (1.0 - current / (memLength * memWidth) / iLimDilCem);
        mc::FFVar concOutConcIntAem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        mc::FFVar concOutDilIntAem = concOutDiluteNd * (1.0 - current / (memLength * memWidth) / iLimDilAem);

        mc::FFVar avgConcConcIntCem = 0.5 * (concOutConcIntCem * scaleFacConc + concInConcIntCem);
        mc::FFVar avgConcDilIntCem = 0.5 * (concOutDilIntCem * scaleFacConc + concInDilIntCem);
        mc::FFVar avgConcConcIntAem = 0.5 * (concOutConcIntAem * scaleFacConc + concInConcIntAem);
        mc::FFVar avgConcDilIntAem = 0.5 * (concOutDilIntAem * scaleFacConc + concInDilIntAem);

        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (current / (memLength * memWidth * faraday));
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConcIntAem - avgConcDilIntAem);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConcIntCem - avgConcDilIntCem);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntAem - osmoticCoeff * avgConcDilIntAem);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntCem - osmoticCoeff * avgConcDilIntCem);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH2O / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        std::vector<mc::FFVar> constraints;
        constraints.reserve(64);

        auto addEqualityConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
            constraints.push_back(-expr);
        };
        auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
        };

        addEqualityConstraint(molNStream1 - molNFeed - molNStream6);
        addEqualityConstraint(molWStream1 - molWFeed - molWStream6);
        addEqualityConstraint(molNStream4 - molNStream3 - molNStream8);
        addEqualityConstraint(molWStream4 - molWStream3 - molWStream8);

        addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1000.0) * molWStream2);
        mc::FFVar purityEd3Expr = molNStream7 - 30000.0 * (molweightH20gMol / molweightN / 1000.0) * molWStream7;
        addLessEqualConstraint(-purityEd3Expr);

        addEqualityConstraint(resConcentrate * conducConcentrate * avgConcConc - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * avgConcDil - thicknessDilute);

        addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);

        mc::FFVar resOhmic = resConcentrate + resDilute + resCem + resAem;
        mc::FFVar voltNonOhm = voltNonOhmicCem + voltNonOhmicAem;

        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInConc) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInDil) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint(voltCellPair - voltNonOhm - resOhmic * (current / (memLength * memWidth)));

        addLessEqualConstraint(current / (memLength * memWidth * iLimDilCem) - 1.0);
        addLessEqualConstraint(current / (memLength * memWidth * iLimDilAem) - 1.0);

        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        mc::FFVar voltBlank = resBlank * current / (memLength * memWidth);
        mc::FFVar voltTotal = voltBlank + numCells * voltCellPair;
        addEqualityConstraint(opex1 - (costElec * (24.0 * 1e-3) * voltTotal * current * daysOperation));

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        mc::FFVar fConc = 1400.0 / reConc;
        mc::FFVar fDil = 104.5 / pow(reDil, 0.37);
        mc::FFVar presDrop1 = densityManure * fConc * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate);
        mc::FFVar presDrop2 = densityManure * fDil * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute);
        mc::FFVar pumpingCost = costElec * (24.0 * 1e-3) *
            (presDrop1 * ((flowOutConc + flowInConc) / 2.0) + presDrop2 * ((flowOutDil + flowInDil) / 2.0)) *
            numCells * daysOperation;
        addEqualityConstraint(opex2 - pumpingCost);

        if (unitIdx == 1) {
            mc::FFVar flowStream2 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream2 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream2 - flowStream2 * concStream2);
            addEqualityConstraint(molWStream2 - flowStream2 * (densityH2OGm3 - (concStream2 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream3 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream3 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream3 - flowStream3 * concStream3);
            addEqualityConstraint(molWStream3 - flowStream3 * (densityH2OGm3 - (concStream3 * molweightN) / molweightH20gMol));
        }

        if (unitIdx == 2) {
            mc::FFVar flowStream5 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream5 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream5 - flowStream5 * concStream5);
            addEqualityConstraint(molWStream5 - flowStream5 * (densityH2OGm3 - (concStream5 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream6 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream6 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream6 - flowStream6 * concStream6);
            addEqualityConstraint(molWStream6 - flowStream6 * (densityH2OGm3 - (concStream6 * molweightN) / molweightH20gMol));
        }

        if (unitIdx == 3) {
            mc::FFVar flowStream7 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream7 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream7 - flowStream7 * concStream7);
            addEqualityConstraint(molWStream7 - flowStream7 * (densityH2OGm3 - (concStream7 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream8 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream8 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream8 - flowStream8 * concStream8);
            addEqualityConstraint(molWStream8 - flowStream8 * (densityH2OGm3 - (concStream8 * molweightN) / molweightH20gMol));
        }

        mc::FFVar objective = this->probability * (capex + opex1 + opex2 + costTotal * 0.0);

        std::vector<mc::FFVar> functions;
        functions.reserve(1 + constraints.size());
        functions.push_back(objective);
        functions.insert(functions.end(), constraints.begin(), constraints.end());
        this->F[scenarioName] = functions;
    }
}

void EDUnits::buildFullModelDAG() {
    const double numCells = 500.0;
    const double daysOperation = 7.0 * 365.0;
    const double rg = 8.314;
    const double temp = 273.15 + 25.0;
    const double faraday = 96485.3321;
    const double molweightN = 14.0067;
    const double molweightNh4Cl = 53.491;
    const double thicknessCem = 0.00011;
    const double thicknessAem = 0.000105;
    const double permSelCem = 0.9;
    const double permSelAem = 0.9;
    const double transCem = 0.95;
    const double transAem = 0.95;
    const double transIonConc = 0.491;
    const double transIonDil = 0.509;
    const double saltDiffCem = 1e-12;
    const double saltDiffAem = 1e-12;
    const double saltDiffConc = 1.84e-9;
    const double saltDiffDil = 1.84e-9;
    const double waterPermCem = 7.79e-11 / 3600.0;
    const double waterPermAem = 6.29e-11 / 3600.0;
    const double vantHoffNumber = 2.0;
    const double molweightH2O = 0.018;
    const double molweightH20gMol = 18.01528;
    const double densityH2O = 1000.0;
    const double densityH2OGm3 = 1000.0 * 1000.0;
    const double densityManure = 1000.0;
    const double viscosityManure = 1e-3;
    const double waterTransNumber = 6.0 + 8.0 + 2.0;
    const double osmoticCoeff = 1.0;
    const double resBlank = 2.5e-5;
    const double resCem = 2.5e-4;
    const double resAem = 1.8e-4;
    const double conducConcentrate = 149.6e-4;
    const double conducDilute = 149.6e-4;
    const double costElec = 0.08147;

    const double flowFeed = 192555.0 / (1000.0 * 24.0 * 60.0 * 60.0);
    const double concFeed = 1366.0;
    const double molNFeed = flowFeed * concFeed / molweightN;
    const double molWFeed = flowFeed * (densityH2OGm3 - concFeed) / molweightH20gMol;

    const double scaleFacFlowByUnit[4] = {
        0.0,
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0),
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0),
        192555.0 / (1000.0 * 24.0 * 60.0 * 60.0)
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        1366.0 / molweightN,
        1366.0 / molweightN,
        1366.0 / molweightN
    };

    const int firstStageCount = static_cast<int>(this->first_stage_IX.size());
    const int secondStageCount = static_cast<int>(this->second_stage_IX.size());
    const int totalSecondStageCount = secondStageCount * static_cast<int>(this->scenario_names.size());

    this->X[ScenarioNames::SCENARIO1].resize(firstStageCount + totalSecondStageCount);
    for (int i = 0; i < firstStageCount + totalSecondStageCount; ++i) {
        this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);
    }

    auto& vars = this->X[ScenarioNames::SCENARIO1];

    mc::FFVar& molNStream1 = vars[0];
    mc::FFVar& molWStream1 = vars[1];
    mc::FFVar& molNStream2 = vars[2];
    mc::FFVar& molWStream2 = vars[3];
    mc::FFVar& molNStream3 = vars[4];
    mc::FFVar& molWStream3 = vars[5];
    mc::FFVar& molNStream4 = vars[6];
    mc::FFVar& molWStream4 = vars[7];
    mc::FFVar& molNStream5 = vars[8];
    mc::FFVar& molWStream5 = vars[9];
    mc::FFVar& molNStream6 = vars[10];
    mc::FFVar& molWStream6 = vars[11];
    mc::FFVar& molNStream7 = vars[12];
    mc::FFVar& molWStream7 = vars[13];
    mc::FFVar& molNStream8 = vars[14];
    mc::FFVar& molWStream8 = vars[15];

    std::vector<mc::FFVar> constraints;
    constraints.reserve(256);

    auto addEqualityConstraint = [&](const mc::FFVar& expr) {
        constraints.push_back(expr);
        constraints.push_back(-expr);
    };
    auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
        constraints.push_back(expr);
    };

    addEqualityConstraint(molNStream1 - molNFeed - molNStream6);
    addEqualityConstraint(molWStream1 - molWFeed - molWStream6);
    addEqualityConstraint(molNStream4 - molNStream3 - molNStream8);
    addEqualityConstraint(molWStream4 - molWStream3 - molWStream8);

    addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1000.0) * molWStream2);
    mc::FFVar purityEd3Expr = molNStream7 - 30000.0 * (molweightH20gMol / molweightN / 1000.0) * molWStream7;
    addLessEqualConstraint(-purityEd3Expr);

    mc::FFVar objective = 0.0;

    for (size_t scenarioIndex = 0; scenarioIndex < this->scenario_names.size(); ++scenarioIndex) {
        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const int secondStageStart = firstStageCount + static_cast<int>(scenarioIndex) * secondStageCount;

        mc::FFVar& current = vars[secondStageStart + 0];
        mc::FFVar& flowSplit = vars[secondStageStart + 1];
        mc::FFVar& memLength = vars[secondStageStart + 2];
        mc::FFVar& memWidth = vars[secondStageStart + 3];
        mc::FFVar& thicknessConcentrate = vars[secondStageStart + 4];
        mc::FFVar& thicknessDilute = vars[secondStageStart + 5];
        mc::FFVar& flowOutConcentrateNd = vars[secondStageStart + 6];
        mc::FFVar& flowOutDiluteNd = vars[secondStageStart + 7];
        mc::FFVar& concOutConcentrateNd = vars[secondStageStart + 8];
        mc::FFVar& concOutDiluteNd = vars[secondStageStart + 9];
        mc::FFVar& voltCellPair = vars[secondStageStart + 10];
        mc::FFVar& resConcentrate = vars[secondStageStart + 11];
        mc::FFVar& resDilute = vars[secondStageStart + 12];
        mc::FFVar& voltNonOhmicCem = vars[secondStageStart + 13];
        mc::FFVar& voltNonOhmicAem = vars[secondStageStart + 14];
        mc::FFVar& capex = vars[secondStageStart + 15];
        mc::FFVar& opex1 = vars[secondStageStart + 16];
        mc::FFVar& opex2 = vars[secondStageStart + 17];
        mc::FFVar& costTotal = vars[secondStageStart + 18];

        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar massStream1 = molNStream1 * molweightN + molWStream1 * molweightH20gMol;
        mc::FFVar flowStream1 = massStream1 / densityH2OGm3;
        mc::FFVar concStream1 = molNStream1 / flowStream1;

        mc::FFVar massStream4 = molNStream4 * molweightN + molWStream4 * molweightH20gMol;
        mc::FFVar flowStream4 = massStream4 / densityH2OGm3;
        mc::FFVar concStream4 = molNStream4 / flowStream4;

        mc::FFVar massStream5 = molNStream5 * molweightN + molWStream5 * molweightH20gMol;
        mc::FFVar flowStream5In = massStream5 / densityH2OGm3;
        mc::FFVar concStream5In = molNStream5 / flowStream5In;

        mc::FFVar flowInEd = flowStream1;
        mc::FFVar concInEd = concStream1;
        if (unitIdx == 2) {
            flowInEd = flowStream4;
            concInEd = concStream4;
        } else if (unitIdx == 3) {
            flowInEd = flowStream5In;
            concInEd = concStream5In;
        }

        mc::FFVar flowInConc = flowSplit * flowInEd / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInEd / numCells;
        mc::FFVar concInConc = concInEd;
        mc::FFVar concInDil = concInEd;

        mc::FFVar uConc = 0.5 * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * thicknessConcentrate);
        mc::FFVar uDil = 0.5 * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * thicknessDilute);

        mc::FFVar reConc = densityManure * uConc * (2.0 * thicknessConcentrate) / viscosityManure;
        mc::FFVar reDil = densityManure * uDil * (2.0 * thicknessDilute) / viscosityManure;

        mc::FFVar scConc = viscosityManure / (densityManure * saltDiffConc);
        mc::FFVar scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConcCem = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDilCem = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);
        mc::FFVar shConcAem = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDilAem = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInConc);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInDil);

        mc::FFVar iLimConcCem = shConcCem * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transCem - transIonConc));
        mc::FFVar iLimDilCem = shDilCem * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transCem - transIonConc));
        mc::FFVar iLimConcAem = shConcAem * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transAem - transIonDil));
        mc::FFVar iLimDilAem = shDilAem * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transAem - transIonDil));

        mc::FFVar concInConcIntCem = concInConc * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concInDilIntCem = concInDil * (1.0 - current / (memLength * memWidth) / iLimDilCem);
        mc::FFVar concInConcIntAem = concInConc * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        mc::FFVar concInDilIntAem = concInDil * (1.0 - current / (memLength * memWidth) / iLimDilAem);

        mc::FFVar concOutConcIntCem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concOutDilIntCem = concOutDiluteNd * (1.0 - current / (memLength * memWidth) / iLimDilCem);
        mc::FFVar concOutConcIntAem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        mc::FFVar concOutDilIntAem = concOutDiluteNd * (1.0 - current / (memLength * memWidth) / iLimDilAem);

        mc::FFVar avgConcConcIntCem = 0.5 * (concOutConcIntCem * scaleFacConc + concInConcIntCem);
        mc::FFVar avgConcDilIntCem = 0.5 * (concOutDilIntCem * scaleFacConc + concInDilIntCem);
        mc::FFVar avgConcConcIntAem = 0.5 * (concOutConcIntAem * scaleFacConc + concInConcIntAem);
        mc::FFVar avgConcDilIntAem = 0.5 * (concOutDilIntAem * scaleFacConc + concInDilIntAem);

        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (current / (memLength * memWidth * faraday));
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConcIntAem - avgConcDilIntAem);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConcIntCem - avgConcDilIntCem);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntAem - osmoticCoeff * avgConcDilIntAem);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntCem - osmoticCoeff * avgConcDilIntCem);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH2O / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        addEqualityConstraint(resConcentrate * conducConcentrate * avgConcConc - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * avgConcDil - thicknessDilute);

        addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);

        mc::FFVar resOhmic = resConcentrate + resDilute + resCem + resAem;
        mc::FFVar voltNonOhm = voltNonOhmicCem + voltNonOhmicAem;

        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInConc) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInDil) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint(voltCellPair - voltNonOhm - resOhmic * (current / (memLength * memWidth)));

        addLessEqualConstraint(current / (memLength * memWidth * iLimDilCem) - 1.0);
        addLessEqualConstraint(current / (memLength * memWidth * iLimDilAem) - 1.0);

        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        mc::FFVar voltBlank = resBlank * current / (memLength * memWidth);
        mc::FFVar voltTotal = voltBlank + numCells * voltCellPair;
        addEqualityConstraint(opex1 - (costElec * (24.0 * 1e-3) * voltTotal * current * daysOperation));

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        mc::FFVar fConc = 1400.0 / reConc;
        mc::FFVar fDil = 104.5 / pow(reDil, 0.37);
        mc::FFVar presDrop1 = densityManure * fConc * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate);
        mc::FFVar presDrop2 = densityManure * fDil * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute);
        mc::FFVar pumpingCost = costElec * (24.0 * 1e-3) *
            (presDrop1 * ((flowOutConc + flowInConc) / 2.0) + presDrop2 * ((flowOutDil + flowInDil) / 2.0)) *
            numCells * daysOperation;
        addEqualityConstraint(opex2 - pumpingCost);

        if (unitIdx == 1) {
            mc::FFVar flowStream2 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream2 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream2 - flowStream2 * concStream2);
            addEqualityConstraint(molWStream2 - flowStream2 * (densityH2OGm3 - (concStream2 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream3 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream3 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream3 - flowStream3 * concStream3);
            addEqualityConstraint(molWStream3 - flowStream3 * (densityH2OGm3 - (concStream3 * molweightN) / molweightH20gMol));
        }

        if (unitIdx == 2) {
            mc::FFVar flowStream5 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream5 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream5 - flowStream5 * concStream5);
            addEqualityConstraint(molWStream5 - flowStream5 * (densityH2OGm3 - (concStream5 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream6 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream6 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream6 - flowStream6 * concStream6);
            addEqualityConstraint(molWStream6 - flowStream6 * (densityH2OGm3 - (concStream6 * molweightN) / molweightH20gMol));
        }

        if (unitIdx == 3) {
            mc::FFVar flowStream7 = flowOutConcentrateNd * scaleFacFlow * numCells;
            mc::FFVar concStream7 = concOutConcentrateNd * scaleFacConc;
            addEqualityConstraint(molNStream7 - flowStream7 * concStream7);
            addEqualityConstraint(molWStream7 - flowStream7 * (densityH2OGm3 - (concStream7 * molweightN) / molweightH20gMol));

            mc::FFVar flowStream8 = flowOutDiluteNd * scaleFacFlow * numCells;
            mc::FFVar concStream8 = concOutDiluteNd * scaleFacConc;
            addEqualityConstraint(molNStream8 - flowStream8 * concStream8);
            addEqualityConstraint(molWStream8 - flowStream8 * (densityH2OGm3 - (concStream8 * molweightN) / molweightH20gMol));
        }

        objective += this->probability * (capex + opex1 + opex2 + costTotal * 0.0);
    }

    std::vector<mc::FFVar> functions;
    functions.reserve(1 + constraints.size());
    functions.push_back(objective);
    functions.insert(functions.end(), constraints.begin(), constraints.end());
    this->F[ScenarioNames::SCENARIO1] = functions;
    this->full_model_built = true;
}

Ipopt::SmartPtr<STModel> EDUnits::clone() {
    Ipopt::SmartPtr<EDUnits> copy = new EDUnits();

    copy->scenario_name = this->scenario_name;
    copy->first_stage_IX = this->first_stage_IX;
    copy->second_stage_IX = this->second_stage_IX;
    copy->scenario_names = this->scenario_names;
    copy->probability = this->probability;
    copy->clearDAG();
    if (this->full_model_built) {
        copy->buildFullModelDAG();
    } else {
        copy->buildDAG();
    }
    return copy;
}