#include "EDUnits.h"
#include <vector>
#include <cmath>


// NEEDED: no division, simplier formation, axu var
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
        mc::Interval(1e-5, 11), // molNStream1
        mc::Interval(1e-5, 500), // molWStream1
        mc::Interval(1e-5, 12), // molNStream2
        mc::Interval(1e-5, 500), // molWStream2
        mc::Interval(1e-5, 13), // molNStream3
        mc::Interval(1e-5, 500), // molWStream3
        mc::Interval(1e-5, 14), // molNStream4
        mc::Interval(1e-5, 500), // molWStream4
        mc::Interval(1e-5, 15), // molNStream5
        mc::Interval(1e-5, 500), // molWStream5
        mc::Interval(1e-5, 16), // molNStream6
        mc::Interval(1e-5, 500), // molWStream6
        mc::Interval(1e-5, 17), // molNStream7
        mc::Interval(1e-5, 500), // molWStream7
        mc::Interval(1e-5, 18), // molNStream8
        mc::Interval(1e-5, 500) // molWStream8
    };
    // this->first_stage_IX = {
    //     mc::Interval(0.222018, 0.222018), // molNStream1
    //     mc::Interval(123.686, 123.686),   // molWStream1
    //     mc::Interval(0.06148, 0.06148),   // molNStream2
    //     mc::Interval(119.5, 119.5),       // molWStream2
    //     mc::Interval(0.160538, 0.160538), // molNStream3
    //     mc::Interval(4.18544, 4.18544),   // molWStream3
    //     mc::Interval(14, 14),             // molNStream4
    //     mc::Interval(365.954, 365.954),   // molWStream4
    //     mc::Interval(13.9953, 13.9953),   // molNStream5
    //     mc::Interval(365.808, 365.808),   // molWStream5
    //     mc::Interval(0.00466994, 0.00466994), // molNStream6
    //     mc::Interval(0.145909, 0.145909), // molWStream6
    //     mc::Interval(0.155868, 0.155868), // molNStream7
    //     mc::Interval(4.03953, 4.03953),   // molWStream7
    //     mc::Interval(13.8395, 13.8395),   // molNStream8
    //     mc::Interval(361.769, 361.769)    // molWStream8
    // };

    // this->second_stage_IX = {
    //     mc::Interval(1, 200),           // I 16
    //     mc::Interval(0, 1),           // flowSplit 17
    //     mc::Interval(0.01, 10),          // memLength 18
    //     mc::Interval(0.01, 10),         // memWidth 19
    //     mc::Interval(0.001, 0.01),     // thicknessConcentrate 20
    //     mc::Interval(0.001, 0.01),     // thicknessDilute 21
    //     mc::Interval(1e-3, 100),         // flowOutConcentrateND 22
    //     mc::Interval(1e-3, 10),         // flowOutDiluteND 23
    //     mc::Interval(1e-3, 100),          // concOutConcentrateND 24
    //     mc::Interval(1e-3, 100),          // concOutDiluteND 25
    //     mc::Interval(0, 5),          // voltCellPair 26
    //     mc::Interval(1e-6, 1e-2),       // resConcentrate 27
    //     mc::Interval(1e-6, 1e-2),       // resDilute 28
    //     mc::Interval(1e-4, 0.25),       // voltNonOhmicCEM 29
    //     mc::Interval(1e-4, 0.25),       // voltNonOhmicAEM 30
    //     mc::Interval(0, 1e6),           // capex
    //     mc::Interval(0, 1e6),           // opex1
    //     mc::Interval(0, 1e6),           // opex2
    //     mc::Interval(0, 1e6),            // costTotal
    //     mc::Interval(1e-3, 1e6),            // auxVar
    //     mc::Interval(1e-3, 1e6)            // auxVar2


    // };  
    this->second_stage_IX = {
        mc::Interval(1, 1+1),                    // I 16
        mc::Interval(0.999486-0.1, 0.999486),      // flowSplit 17
        mc::Interval(0.511088, 0.511088+0.1),      // memLength 18
        mc::Interval(0.00999999, 0.00999999+0.1),  // memWidth 19
        mc::Interval(0.00099999, 0.00099999+0.1),  // thicknessConcentrate 20
        mc::Interval(0.00099999, 0.00099999+0.1),  // thicknessDilute 21
        mc::Interval(3.39017, 3.39017+0.1),        // flowOutConcentrateND 22
        mc::Interval(0.00099999, 0.00099999+0.1),  // flowOutDiluteND 23
        mc::Interval(2.52589, 2.52589+0.1),        // concOutConcentrateND 24
        mc::Interval(0.0719748, 0.0719748+0.1),    // concOutDiluteND 25
        mc::Interval(0.206921, 0.206921+0.1),      // voltCellPair 26
        mc::Interval(6.6719e-05, 6.6719e-05+0.1),  // resConcentrate 27
        mc::Interval(0.000129758, 0.000129758+0.1),// resDilute 28
        mc::Interval(0.0431747, 0.0431747+0.1),    // voltNonOhmicCEM 29
        mc::Interval(0.0411688, 0.0411688+0.1),    // voltNonOhmicAEM 30
        mc::Interval(543.797, 543.797+0.1),        // capex
        mc::Interval(516.886, 516.886+0.1),        // opex1
        mc::Interval(0.00745432, 0.00745432+0.1),  // opex2
        mc::Interval(500000, 500000+0.1),          // costTotal
        mc::Interval(155.954, 155.954+0.1),        // auxVar
        mc::Interval(170.04, 170.04+0.1)           // auxVar2
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
    const double molweightH20 = 0.018;
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
        5.388e-06,
        5.388e-06,
        5.388e-06
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        3.967e+02,
        3.967e+02,
        3.967e+02
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
        mc::FFVar& avgConcDilIntCem = vars[secondStageStart + 19];
        mc::FFVar& avgConcDilIntAem = vars[secondStageStart + 20];


        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream1*densityH2OGm3 / (molNStream1*molweightN + molWStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream4*densityH2OGm3 / (molNStream4*molweightN + molWStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream5*densityH2OGm3 / (molNStream5*molweightN + molWStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;

        mc::FFVar reConc = densityManure * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * viscosityManure);
        mc::FFVar reDil = densityManure * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * viscosityManure);

        const double scConc = viscosityManure / (densityManure * saltDiffConc);
        const double scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConc = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDil = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInED);

        // avgConc*Int* directly, avgConcConc/avgConcDil cancelled out of the correction term
        mc::FFVar avgConcConcIntCem = avgConcConc +
            (current / (memLength * memWidth)) * (2.0 * thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc);

        mc::FFVar avgConcConcIntAem = avgConcConc +
            (current / (memLength * memWidth)) * (2.0 * thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc);




        // mc::FFVar avgConcDilIntCem = avgConcDil -
        //     (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil);

        // mc::FFVar avgConcDilIntAem = avgConcDil -
        //     (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil);
        
        
        
        
        std::vector<mc::FFVar> constraints;
        auto addEqualityConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
            constraints.push_back(-expr);
        };
        auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
        };

        addEqualityConstraint( avgConcDilIntCem -avgConcDil +
            (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil));

        addEqualityConstraint( avgConcDilIntAem - avgConcDil +
            (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil));
        
        
        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (current / (memLength * memWidth * faraday));
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConcIntAem - avgConcDilIntAem);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConcIntCem - avgConcDilIntCem);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * osmoticCoeff * (avgConcConcIntAem - avgConcDilIntAem);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * osmoticCoeff * (avgConcConcIntCem - avgConcDilIntCem);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;



        addEqualityConstraint(molNStream1 - molNFeed - molNStream6);
        addEqualityConstraint(molWStream1 - molWFeed - molWStream6);
        addEqualityConstraint(molNStream4 - molNStream3 - molNStream8);
        addEqualityConstraint(molWStream4 - molWStream3 - molWStream8);

        addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1e6) * molWStream2);
        addLessEqualConstraint(30000.0 * (molweightH20gMol / molweightN / 1e6) * molWStream7 - molNStream7);

        addLessEqualConstraint(molNStream1 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream1 * molweightH20gMol);
        addLessEqualConstraint(molNStream2 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream2 * molweightH20gMol);
        addLessEqualConstraint(molNStream3 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream3 * molweightH20gMol);
        addLessEqualConstraint(molNStream4 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream4 * molweightH20gMol);
        addLessEqualConstraint(molNStream5 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream5 * molweightH20gMol);
        addLessEqualConstraint(molNStream6 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream6 * molweightH20gMol);
        addLessEqualConstraint(molNStream7 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream7 * molweightH20gMol);
        addLessEqualConstraint(molNStream8 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream8 * molweightH20gMol);

        addEqualityConstraint(resConcentrate * conducConcentrate * avgConcConc - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * avgConcDil - thicknessDilute);


        //addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        //addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);


        addEqualityConstraint(voltNonOhmicCem * (faraday / (permSelCem * rg * temp)) - log(avgConcConcIntCem)+log(avgConcDilIntCem));



        addEqualityConstraint(voltNonOhmicAem * (faraday / (permSelAem * rg * temp)) - log(avgConcConcIntAem) +log(avgConcDilIntAem));


        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));
        

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));



        // REDUNDANT CONSTRAINT: flow balance is already enforced by molNStream and molWStream constraints, but this may help the solver
        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow)+flowOutDiluteNd - (flowInDil / scaleFacFlow));
        //Redundant constraint: ion flux is already enforced by flowOutConc/Dil*concOutConc/Dil constraints, but may help solver
        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc)+(flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc));


            
        addEqualityConstraint(voltCellPair - (voltNonOhmicCem + voltNonOhmicAem) -
            (resConcentrate + resDilute + resCem + resAem) * (current / (memLength * memWidth)));


        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * current / (memLength * memWidth) + numCells * voltCellPair) * current * daysOperation);

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
           densityManure * (1400.0 / reConc) * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate) * ((flowOutConc + flowInConc) / 2.0) +
           densityManure * (104.5 / pow(reDil, 0.37)) * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute) * ((flowOutDil + flowInDil) / 2.0)));

        if (unitIdx == 1) {
            addEqualityConstraint(molNStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 2) {
            addEqualityConstraint(molNStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 3) {
            addEqualityConstraint(molNStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        mc::FFVar objective = this->probability * (capex + opex1 + opex2);

        std::vector<mc::FFVar> functions;

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
    const double molweightH20 = 0.018;
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
        5.388e-06,
        5.388e-06,
        5.388e-06
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        3.967e+02,
        3.967e+02,
        3.967e+02
    };

    const int firstStageCount = static_cast<int>(this->first_stage_IX.size());
    // second_stage_IX is already the full combined size after convertToCentralizedModel()
    const int totalSecondStageCount = static_cast<int>(this->second_stage_IX.size());
    const int secondStageCount = totalSecondStageCount / static_cast<int>(this->scenario_names.size());

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

    addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1e6) * molWStream2);
    addLessEqualConstraint(30000.0 * (molweightH20gMol / molweightN / 1e6) * molWStream7 - molNStream7);

    addLessEqualConstraint(molNStream1 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream1 * molweightH20gMol);
    addLessEqualConstraint(molNStream2 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream2 * molweightH20gMol);
    addLessEqualConstraint(molNStream3 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream3 * molweightH20gMol);
    addLessEqualConstraint(molNStream4 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream4 * molweightH20gMol);
    addLessEqualConstraint(molNStream5 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream5 * molweightH20gMol);
    addLessEqualConstraint(molNStream6 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream6 * molweightH20gMol);
    addLessEqualConstraint(molNStream7 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream7 * molweightH20gMol);
    addLessEqualConstraint(molNStream8 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream8 * molweightH20gMol);

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
        mc::FFVar& avgConcDilIntAem = vars[secondStageStart + 19];
        mc::FFVar& avgConcDilIntCem = vars[secondStageStart + 20];

 

        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED =  densityH2OGm3 / (molweightN + molWStream1/molNStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = densityH2OGm3 / (molweightN + molWStream4/molNStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = densityH2OGm3 / (molweightN + molWStream5/molNStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;

        mc::FFVar reConc = densityManure * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * viscosityManure);
        mc::FFVar reDil = densityManure * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * viscosityManure);

        const double scConc = viscosityManure / (densityManure * saltDiffConc);
        const double scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConc = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDil = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInED);

        mc::FFVar iLimConcCem = shConc * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transCem - transIonConc));
        mc::FFVar iLimDilCem = shDil * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transCem - transIonConc));
        mc::FFVar iLimConcAem = shConc * avgConcConc * faraday * saltDiffConc / (2.0 * thicknessConcentrate * (transAem - transIonDil));
        mc::FFVar iLimDilAem = shDil * avgConcDil * faraday * saltDiffDil / (2.0 * thicknessDilute * (transAem - transIonDil));

        mc::FFVar concInConcIntCem = concInED * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concInConcIntAem = concInED * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        mc::FFVar concOutConcIntCem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcCem);
        mc::FFVar concOutConcIntAem = concOutConcentrateNd * (1.0 + current / (memLength * memWidth) / iLimConcAem);
        



        mc::FFVar concInDilIntCem  = concInED        * (1.0 - current / (memLength * memWidth * iLimDilCem));
        mc::FFVar concOutDilIntCem = concOutDiluteNd * (1.0 - current / (memLength * memWidth * iLimDilCem));
        mc::FFVar concInDilIntAem  = concInED        * (1.0 - current / (memLength * memWidth * iLimDilAem));
        mc::FFVar concOutDilIntAem = concOutDiluteNd * (1.0 - current / (memLength * memWidth * iLimDilAem));


        mc::FFVar avgConcConcIntCem = 0.5 * (concOutConcIntCem * scaleFacConc + concInConcIntCem);
        mc::FFVar avgConcConcIntAem = 0.5 * (concOutConcIntAem * scaleFacConc + concInConcIntAem);


        addEqualityConstraint(avgConcDilIntCem - 0.5 * (concOutDilIntCem * scaleFacConc + concInDilIntCem));
        addEqualityConstraint(avgConcDilIntAem - 0.5 * (concOutDilIntAem * scaleFacConc + concInDilIntAem));


        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (current / (memLength * memWidth * faraday));
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConcIntAem - avgConcDilIntAem);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConcIntCem - avgConcDilIntCem);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntAem - osmoticCoeff * avgConcDilIntAem);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConcIntCem - osmoticCoeff * avgConcDilIntCem);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        addEqualityConstraint(resConcentrate * conducConcentrate * avgConcConc - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * avgConcDil - thicknessDilute);


        //addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        //addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);
        addEqualityConstraint(
            voltNonOhmicCem * (faraday / (permSelCem * rg * temp)) 
            - log(avgConcConcIntCem / avgConcDilIntCem)
        );

        addEqualityConstraint(
            voltNonOhmicAem * (faraday / (permSelAem * rg * temp)) 
            - log(avgConcConcIntAem / avgConcDilIntAem)
        );
        
       // fix : 1/x vs log(x) issues
        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));


        // REDUNDANT CONSTRAINT: flow balance is already enforced by molNStream and molWStream constraints, but this may help the solver
        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow)+flowOutDiluteNd - (flowInDil / scaleFacFlow));
        //Redundant constraint: ion flux is already enforced by flowOutConc/Dil*concOutConc/Dil constraints, but may help solver
        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc)+(flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc));



        addEqualityConstraint(voltCellPair - (voltNonOhmicCem + voltNonOhmicAem) -
            (resConcentrate + resDilute + resCem + resAem) * (current / (memLength * memWidth)));

        addLessEqualConstraint(current / (memLength * memWidth * iLimDilCem) - 1.0);
        addLessEqualConstraint(current / (memLength * memWidth * iLimDilAem) - 1.0);

        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * current / (memLength * memWidth) + numCells * voltCellPair) * current * daysOperation);

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
            densityManure * (1400.0 / reConc) * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate) * ((flowOutConc + flowInConc) / 2.0) +
            densityManure * (104.5 / pow(reDil, 0.37)) * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute) * ((flowOutDil + flowInDil) / 2.0)));

        if (unitIdx == 1) {
            addEqualityConstraint(molNStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 2) {
            addEqualityConstraint(molNStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 3) {
            addEqualityConstraint(molNStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        objective += this->probability * (capex + opex1 + opex2 + costTotal * 0.0);
    }

    std::vector<mc::FFVar> functions;

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