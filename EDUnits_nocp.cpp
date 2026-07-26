#include "EDUnits_nocp.h"

#include <vector>
#include <cmath>

EDUnits_nocp::EDUnits_nocp(BranchingStrategy branchingStrategy) : STModel() {
    this->branching_strategy = branchingStrategy;
    this->scenario_names = {
        ScenarioNames::SCENARIO1,
        ScenarioNames::SCENARIO2,
        ScenarioNames::SCENARIO3
    };
    this->scenario_name = ScenarioNames::SCENARIO1;
    this->probability = 1.0;

    this->first_stage_IX = {
        mc::Interval(0.22127814806092042, 0.22127814806092042), // molNStream1
        mc::Interval(123.68614069261695, 123.68614069261695), // molWStream1
        mc::Interval(0.06147995475275122, 0.06147995475275122), // molNStream2
        mc::Interval(119.50012464909791, 119.50012464909791), // molWStream2
        mc::Interval(0.15979819330816922, 0.15979819330816922), // molNStream3
        mc::Interval(4.186016043519042, 4.186016043519042), // molWStream3
        mc::Interval(0.19509784106870087, 0.19509784106870087), // molNStream4
        mc::Interval(5.156888893299397, 5.156888893299397), // molWStream4
        mc::Interval(0.19116783413547886, 0.19116783413547886), // molNStream5
        mc::Interval(5.010404710323144, 5.010404710323144), // molWStream5
        mc::Interval(0.0039300069332220125, 0.0039300069332220125), // molNStream6
        mc::Interval(0.14648418297625343, 0.14648418297625343), // molWStream6
        mc::Interval(0.1558681863749472, 0.1558681863749472), // molNStream7
        mc::Interval(4.03953186054279, 4.03953186054279), // molWStream7
        mc::Interval(0.03529964776053168, 0.03529964776053168), // molNStream8
        mc::Interval(0.9708728497803544, 0.9708728497803544) // molWStream8
    };
    // this->first_stage_IX = {
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500),
    //     mc::Interval(1e-5, 10),
    //     mc::Interval(1e-5, 500)
    // };

    this->second_stage_IX = {
        mc::Interval(200,200),           // current 16
        mc::Interval(0, 1),           // flowSplit 17
        mc::Interval(10, 10),          // memLength 18
        mc::Interval(0.01, 10),         // memWidth 19
        mc::Interval(0.001, 0.01),     // thicknessConcentrate 20
        mc::Interval(0.001, 0.01),     // thicknessDilute 21
        mc::Interval(1e-3, 100),         // flowOutConcentrateND 22
        mc::Interval(1e-3, 10),         // flowOutDiluteND 23
        mc::Interval(1e-3, 100),          // concOutConcentrateND 24
        mc::Interval(1e-3, 100),          // concOutDiluteND 25
        mc::Interval(0, 5),          // voltCellPair 26
        mc::Interval(1e-6, 1e-2),       // resConcentrate 27 
        mc::Interval(1e-6, 1e-2),       // resDilute 28
        mc::Interval(1e-4, 0.25),       // voltNonOhmicCEM 29
        mc::Interval(1e-4, 0.25),       // voltNonOhmicAEM 30
        mc::Interval(0, 1e6),           // capex 31
        mc::Interval(0, 1e6),           // opex1 32
        mc::Interval(0, 1e6),           // opex2 33
        mc::Interval(10, 1e4),           // currentDensity 34


    };
    // this->second_stage_IX = {

    //     mc::Interval (100.5,100.5), // x[16]
    //     mc::Interval (0.956123,0.956123), // x[17]
    //     mc::Interval (5.005,5.005), // x[18]
    //     mc::Interval (0.0550073,0.0550073), // x[19]
    //     mc::Interval (0.00100003,0.00100003), // x[20]
    //     mc::Interval (0.00100002,0.00100002), // x[21]
    //     mc::Interval (0.0344994,0.0344994), // x[22]
    //     mc::Interval (0.001,0.001), // x[23]
    //     mc::Interval (5.18493,5.18493), // x[24]
    //     mc::Interval (3.67734,3.67734), // x[25]
    //     mc::Interval (0.171929,0.171929), // x[26]
    //     mc::Interval (3.26329e-05,3.26329e-05), // x[27]
    //     mc::Interval (3.82108e-05,3.82108e-05), // x[28]
    //     mc::Interval (0.00364867,0.00364867), // x[29]
    //     mc::Interval (0.00364867,0.00364867), // x[30]
    //     mc::Interval (321.946,321.946), // x[31]
    //     mc::Interval (427.181,427.181), // x[32]
    //     mc::Interval (0.0148602,0.0148602), // x[33]
    //     mc::Interval (328.708,328.708), // x[34]
    // };
    }
void EDUnits_nocp::buildDAG() {
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
    const double waterTransNumber = 6.0 + 8.0;
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
        mc::FFVar& currentDensity = vars[secondStageStart + 18];



        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream1 * densityH2OGm3 / (molNStream1 * molweightN + molWStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream4 * densityH2OGm3 / (molNStream4 * molweightN + molWStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream5 * densityH2OGm3 / (molNStream5 * molweightN + molWStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;

        mc::FFVar uConc = 0.5 * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * thicknessConcentrate);
        mc::FFVar uDil = 0.5 * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * thicknessDilute);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInED);

        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * currentDensity / faraday;
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConc - avgConcDil);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConc - avgConcDil);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConc - osmoticCoeff * avgConcDil);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConc - osmoticCoeff * avgConcDil);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        std::vector<mc::FFVar> constraints;
  
        auto addEqualityConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
            constraints.push_back(-expr);
        };
        auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
        };




        // // convex envelope
        // Facet 1 (eq. 3.17) — all lower bounds
        // addLessEqualConstraint(memLength_L*memWidth_L*currentDensity
        // + currentDensity_L*memWidth_L*memLength
        // + currentDensity_L*memLength_L*memWidth
        // - 2*currentDensity_L*memLength_L*memWidth_L-current);

        // // Facet 2 (eq. 3.18) — all upper bounds
        // addLessEqualConstraint(memLength_U*memWidth_U*currentDensity
        // + currentDensity_U*memWidth_U*memLength
        // + currentDensity_U*memLength_U*memWidth
        // - 2*currentDensity_U*memLength_U*memWidth_U-current);

        // // Facet 3 (eq. 3.19)
        // addLessEqualConstraint(memLength_L*memWidth_U*currentDensity
        // + currentDensity_L*memWidth_U*memLength
        // + currentDensity_U*memLength_L*memWidth
        // - currentDensity_U*memLength_L*memWidth_U
        // - currentDensity_L*memLength_L*memWidth_U-current);

        // // Facet 4 (eq. 3.20)
        // addLessEqualConstraint(memLength_U*memWidth_L*currentDensity
        // + currentDensity_U*memWidth_L*memLength
        // + currentDensity_L*memLength_U*memWidth
        // - currentDensity_U*memLength_U*memWidth_L
        // - currentDensity_L*memLength_U*memWidth_L-current);

        // concave envelope
        // addLessEqualConstraint(current-(memLength_L*memWidth_L*currentDensity
        //     + currentDensity_U*memWidth_U*memLength
        //     + currentDensity_U*memLength_L*memWidth
        //     - currentDensity_U*memLength_L*(memWidth_L + memWidth_U)));

        // addLessEqualConstraint(current- (memLength_L*memWidth_U*currentDensity
        //     + currentDensity_U*memWidth_U*memLength
        //     + currentDensity_L*memLength_L*memWidth
        //     - memLength_L*memWidth_U*(currentDensity_U + currentDensity_L)));

        // addLessEqualConstraint(current- (memLength_L*memWidth_L*currentDensity
        //     + currentDensity_U*memWidth_L*memLength
        //     + currentDensity_U*memLength_U*memWidth
        //     - currentDensity_U*memWidth_L*(memLength_U + memLength_L)));

        // addLessEqualConstraint(current- (memLength_U*memWidth_L*currentDensity
        //     + currentDensity_L*memWidth_L*memLength
        //     + currentDensity_U*memLength_U*memWidth
        //     - memLength_U*memWidth_L*(currentDensity_L + currentDensity_U)));

        // addLessEqualConstraint(current- (memLength_U*memWidth_U*currentDensity
        //     + currentDensity_L*memWidth_L*memLength
        //     + currentDensity_L*memLength_U*memWidth
        //     - currentDensity_L*memLength_U*(memWidth_L + memWidth_U)));

        // addLessEqualConstraint(current- (memLength_U*memWidth_U*currentDensity
        //     + currentDensity_L*memWidth_U*memLength
        //     + currentDensity_L*memLength_L*memWidth
        //     - currentDensity_L*memWidth_U*(memLength_L + memLength_U)));



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


        addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDil - avgConcConc);
        addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDil - avgConcConc);


        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * pow(memWidth, 2.0) / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * pow(memWidth, 2.0) / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * pow(memWidth, 2.0) / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * pow(memWidth, 2.0) / scaleFacFlow));

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc) -
            (fluxIonsTotal * pow(memWidth, 2.0) / (scaleFacFlow * scaleFacConc)));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * pow(memWidth, 2.0) / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint(
            voltCellPair - (voltNonOhmicCem + voltNonOhmicAem)
            - (resConcentrate + resDilute + resCem + resAem) * currentDensity
        );

        addEqualityConstraint(capex - (6800.0  + 2.0 * 100.0 * (numCells - 2.0)) * pow(memWidth, 2.0));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * pow(currentDensity, 2.0) + numCells * voltCellPair*currentDensity) * pow(memWidth, 2.0) * daysOperation);

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;

        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
            175.0 * viscosityManure * pow(memWidth, 2.0) * pow(uConc, 2.0) / thicknessConcentrate
            + 78.70934593096298 * pow(memWidth, 1.37)* pow(uDil, 2.0)
                * pow(flowOutDil + flowInDil, 0.63) / thicknessDilute));

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

        mc::FFVar objective = this->probability * (capex + opex1 + opex2 );


        std::vector<mc::FFVar> functions;

        functions.push_back(objective);
        functions.insert(functions.end(), constraints.begin(), constraints.end());
        this->F[scenarioName] = functions;
    }
}

void EDUnits_nocp::buildFullModelDAG() {
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
    const double waterTransNumber = 6.0 + 8.0;
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
        // mc::FFVar& ratioCurrentILimCem   = vars[secondStageStart + 19];  // replaces auxConcDilCem
        // mc::FFVar& ratioCurrentILimAem   = vars[secondStageStart + 20];  // replaces auxConcDilAem

        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream1 * densityH2OGm3 / (molNStream1 * molweightN + molWStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream4 * densityH2OGm3 / (molNStream4 * molweightN + molWStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream5 * densityH2OGm3 / (molNStream5 * molweightN + molWStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;

        mc::FFVar uConc = 0.5 * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * thicknessConcentrate);
        mc::FFVar uDil = 0.5 * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * thicknessDilute);

        mc::FFVar reConc = densityManure * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * viscosityManure);
        mc::FFVar reDil = densityManure * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * viscosityManure);

        mc::FFVar avgConcConc = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED);
        mc::FFVar avgConcDil = 0.5 * (concOutDiluteNd * scaleFacConc + concInED);

        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (current / (memLength * memWidth * faraday));
        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (avgConcConc - avgConcDil);
        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (avgConcConc - avgConcDil);
        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConc - osmoticCoeff * avgConcDil);
        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * (osmoticCoeff * avgConcConc - osmoticCoeff * avgConcDil);
        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        addEqualityConstraint(resConcentrate * conducConcentrate * avgConcConc - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * avgConcDil - thicknessDilute);



        addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDil- avgConcConc);
        addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDil - avgConcConc);



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

        addEqualityConstraint(voltCellPair - (voltNonOhmicCem + voltNonOhmicAem) -
            (resConcentrate + resDilute + resCem + resAem) * (current / (memLength * memWidth)));

        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * current / (memLength * memWidth) + numCells * voltCellPair) * current * daysOperation);

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
            densityManure * (1400.0 / reConc) * memLength * pow(uConc, 2.0) / (4.0 * thicknessConcentrate) * ((flowOutConc + flowInConc) / 2.0) +
            densityManure * (104.5 / pow(reDil, 0.37)) * memLength * pow(uDil, 2.0) / (4.0 * thicknessDilute) * ((flowOutDil + flowInDil) / 2.0)));


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

        objective += this->probability * (capex + opex1 + opex2 );
    }

    std::vector<mc::FFVar> functions;

    functions.push_back(objective);
    functions.insert(functions.end(), constraints.begin(), constraints.end());
    this->F[ScenarioNames::SCENARIO1] = functions;
    this->full_model_built = true;
}

Ipopt::SmartPtr<STModel> EDUnits_nocp::clone() {
    Ipopt::SmartPtr<EDUnits_nocp> copy = new EDUnits_nocp();

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