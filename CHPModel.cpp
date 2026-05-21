#include "CHPModel.h"
#include <cmath>

// Scenario data: { Qdot_dem [MW], P_dem [MW] }
// DEFAULT_DATA from CHP_sizing_SP.py / maingo example
static const std::vector<std::vector<double>> CHP_SCENARIO_DATA = {
    {1.1, 0.5}, {1.1, 2.5}, {1.1, 4.0},
    {1.2, 0.5}, {1.2, 4.0},
    {1.3, 0.5}, {1.3, 2.5}, {1.3, 4.0},
};

// Variable index offsets within a scenario subproblem (Nx=1)
static const int IDX_QDOT_NOM = 0; // first-stage
static const int IDX_QDOT_REL = 1; // second-stage offset 0
static const int IDX_EDOT_IN  = 2; // second-stage offset 1
static const int IDX_P_OUT    = 3; // second-stage offset 2
static const int IDX_P_BUY    = 4; // second-stage offset 3
static const int IDX_P_SELL   = 5; // second-stage offset 4

CHPModel::CHPModel(BranchingStrategy branching_strategy) : STModel() {
    this->branching_strategy = branching_strategy;

    this->scenario_names = {
        ScenarioNames::SCENARIO1, ScenarioNames::SCENARIO2,
        ScenarioNames::SCENARIO3, ScenarioNames::SCENARIO4,
        ScenarioNames::SCENARIO5, ScenarioNames::SCENARIO6,
        ScenarioNames::SCENARIO7, ScenarioNames::SCENARIO8,
    };
    this->scenario_name = ScenarioNames::SCENARIO1;
    this->probability = 1.0 / static_cast<double>(this->scenario_names.size());

    // Store scenario data in perturb_coeffs: [Qdot_dem, P_dem]
    for (int s = 0; s < static_cast<int>(this->scenario_names.size()); ++s) {
        this->perturb_coeffs[this->scenario_names[s]] = CHP_SCENARIO_DATA[s];
    }

    // Problem constants
    const double n_years = 30.0, i_rate = 0.05, c_m = 0.05, c_ref = 1.3e6;
    const double pown = std::pow(1.0 + i_rate, n_years);
    const double af = pown * i_rate / (pown - 1.0);
    this->invest_coeff = 1e-6 * (af + c_m) * c_ref;

    // Bounds
    const double Qdot_nom_min = 1.4, Qdot_nom_max = 2.3;
    const double Qdot_rel_min_bound = 0.5, Qdot_eps = 0.001 * Qdot_nom_max;

    // eff_th(Qdot_nom_max, 1.0) = eta_th_nom(2.3) * eta_th_rel(1.0)
    const double eta_th_max = (0.498 - Qdot_nom_max / 21.17) *
                              (-0.0768 - 0.0199 + 1.0960);
    // eff_el(Qdot_nom_max, 1.0) = eta_el_nom(2.3) * eta_el_rel(1.0)
    const double eta_el_max = (0.372 + Qdot_nom_max / 21.17) *
                              (-0.2611 + 0.6743 + 0.5868);

    const double Edot_in_max = Qdot_nom_max / eta_th_max;
    const double P_out_max   = Edot_in_max * eta_el_max;
    const double P_dem_max   = 4.0; // max P_dem across all scenarios

    // Nx=1: Qdot_nom
    this->first_stage_IX = {
        mc::Interval(Qdot_nom_min, Qdot_nom_max),
    };

    // Ny=5 per scenario: Qdot_rel, Edot_in, P_out, P_buy, P_sell
    this->second_stage_IX = {
        mc::Interval(0.0, 1.0),          // Qdot_rel
        mc::Interval(0.0, Edot_in_max),  // Edot_in
        mc::Interval(0.0, P_out_max),    // P_out
        mc::Interval(0.0, P_dem_max),    // P_buy
        mc::Interval(0.0, P_out_max),    // P_sell
    };
}

void CHPModel::buildDAG() {
    const int Nx = static_cast<int>(this->first_stage_IX.size());  // 1
    const int Ny = static_cast<int>(this->second_stage_IX.size()); // 5
    const int nvars = Nx + Ny;

    // Part-load constraint coefficients (vertex_form expanded to ax²+bx+c)
    // vertex_form(Qdot_rel, a=-1, b=Qdot_rel_min+Qdot_eps, c=-Qdot_rel_min*Qdot_eps)
    const double Qdot_rel_min = 0.5, Qdot_eps = 0.001 * 2.3;
    const double b_part = Qdot_rel_min + Qdot_eps;         // 0.5023
    const double c_part = -Qdot_rel_min * Qdot_eps;        // -0.00115
    const double T_OP_coeff = 1e-6 * 6000.0;               // 1e-6 * T_OP

    for (auto& sn : this->scenario_names) {
        this->X[sn].resize(nvars);
        for (int i = 0; i < nvars; ++i)
            this->X[sn][i].set(&this->DAG[sn]);

        mc::FFVar& Qdot_nom = this->X[sn][IDX_QDOT_NOM];
        mc::FFVar& Qdot_rel = this->X[sn][IDX_QDOT_REL];
        mc::FFVar& Edot_in  = this->X[sn][IDX_EDOT_IN];
        mc::FFVar& P_out    = this->X[sn][IDX_P_OUT];
        mc::FFVar& P_buy    = this->X[sn][IDX_P_BUY];
        mc::FFVar& P_sell   = this->X[sn][IDX_P_SELL];

        double Qdot_dem = this->perturb_coeffs[sn][0];
        double P_dem    = this->perturb_coeffs[sn][1];

        mc::FFVar Qdot_rel_sq = Qdot_rel * Qdot_rel;

        mc::FFVar eta_th = (0.498 - Qdot_nom / 21.17) *
                           (-0.0768 * Qdot_rel_sq - 0.0199 * Qdot_rel + 1.0960);
        mc::FFVar eta_el = (0.372 + Qdot_nom / 21.17) *
                           (-0.2611 * Qdot_rel_sq + 0.6743 * Qdot_rel + 0.5868);

        // Edot_in * eta_th - Qdot_nom * Qdot_rel == 0
        mc::FFVar c_edot = Edot_in * eta_th - Qdot_nom * Qdot_rel;

        // P_out - Edot_in * eta_el == 0
        mc::FFVar c_pout = P_out - Edot_in * eta_el;

        // P_buy - P_sell - P_dem + P_out == 0
        mc::FFVar c_pbal = P_buy - P_sell - P_dem + P_out;

        // Qdot_dem - Qdot_nom * Qdot_rel <= 0
        mc::FFVar c_heat = Qdot_dem - Qdot_nom * Qdot_rel;

        // Part-load: -Qdot_rel^2 + b_part*Qdot_rel + c_part <= 0
        // Feasible when Qdot_rel <= Qdot_eps (off) or Qdot_rel >= Qdot_rel_min (on)
        mc::FFVar c_partload = -Qdot_rel_sq + b_part * Qdot_rel + c_part;

        mc::FFVar obj = this->probability * (
            this->invest_coeff * pow(Qdot_nom, 0.9)
            + T_OP_coeff * (80.0 * Edot_in + 250.0 * P_buy - 100.0 * P_sell)
        );

        this->F[sn] = {obj, c_edot, -c_edot, c_pout, -c_pout, c_pbal, -c_pbal, c_heat, c_partload};
    }
}

void CHPModel::buildFullModelDAG() {
    const int Nx = static_cast<int>(this->first_stage_IX.size());  // 1
    const int Ny = static_cast<int>(this->second_stage_IX.size()); // 5 per scenario
    const int Ns = static_cast<int>(this->scenario_names.size());  // 8

    const double b_part = 0.5 + 0.001 * 2.3;        // 0.5023
    const double c_part = -0.5 * 0.001 * 2.3;        // -0.00115
    const double T_OP_coeff = 1e-6 * 6000.0;

    this->X[ScenarioNames::SCENARIO1].resize(Nx + Ns * Ny);

    for (int i = 0; i < Nx; ++i)
        this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);
    for (int s = 0; s < Ns; ++s) {
        int y0 = Nx + s * Ny;
        for (int i = 0; i < Ny; ++i)
            this->X[ScenarioNames::SCENARIO1][y0 + i].set(&this->DAG[ScenarioNames::SCENARIO1]);
    }

    mc::FFVar objective = 0;

    for (int s = 0; s < Ns; ++s) {
        int y0 = Nx + s * Ny;
        ScenarioNames sn = this->scenario_names[s];

        mc::FFVar& Qdot_nom = this->X[ScenarioNames::SCENARIO1][0];
        mc::FFVar& Qdot_rel = this->X[ScenarioNames::SCENARIO1][y0 + 0];
        mc::FFVar& Edot_in  = this->X[ScenarioNames::SCENARIO1][y0 + 1];
        mc::FFVar& P_out    = this->X[ScenarioNames::SCENARIO1][y0 + 2];
        mc::FFVar& P_buy    = this->X[ScenarioNames::SCENARIO1][y0 + 3];
        mc::FFVar& P_sell   = this->X[ScenarioNames::SCENARIO1][y0 + 4];

        double Qdot_dem = this->perturb_coeffs[sn][0];
        double P_dem    = this->perturb_coeffs[sn][1];

        mc::FFVar Qdot_rel_sq = Qdot_rel * Qdot_rel;

        mc::FFVar eta_th = (0.498 - Qdot_nom / 21.17) *
                           (-0.0768 * Qdot_rel_sq - 0.0199 * Qdot_rel + 1.0960);
        mc::FFVar eta_el = (0.372 + Qdot_nom / 21.17) *
                           (-0.2611 * Qdot_rel_sq + 0.6743 * Qdot_rel + 0.5868);

        mc::FFVar c_edot    = Edot_in * eta_th - Qdot_nom * Qdot_rel;
        mc::FFVar c_pout    = P_out - Edot_in * eta_el;
        mc::FFVar c_pbal    = P_buy - P_sell - P_dem + P_out;
        mc::FFVar c_heat    = Qdot_dem - Qdot_nom * Qdot_rel;
        mc::FFVar c_partload = -Qdot_rel_sq + b_part * Qdot_rel + c_part;

        objective += this->probability * (
            this->invest_coeff * pow(Qdot_nom, 0.9)
            + T_OP_coeff * (80.0 * Edot_in + 250.0 * P_buy - 100.0 * P_sell)
        );

        std::vector<mc::FFVar> scenario_constraints = {
            c_edot, -c_edot, c_pout, -c_pout, c_pbal, -c_pbal, c_heat, c_partload
        };
        this->F[ScenarioNames::SCENARIO1].insert(
            this->F[ScenarioNames::SCENARIO1].end(),
            scenario_constraints.begin(), scenario_constraints.end()
        );
    }

    this->F[ScenarioNames::SCENARIO1].insert(
        this->F[ScenarioNames::SCENARIO1].begin(), objective
    );
    this->full_model_built = true;
}

Ipopt::SmartPtr<STModel> CHPModel::clone() {
    Ipopt::SmartPtr<CHPModel> p = new CHPModel();
    p->scenario_name    = this->scenario_name;
    p->first_stage_IX   = this->first_stage_IX;
    p->second_stage_IX  = this->second_stage_IX;
    p->perturb          = this->perturb;
    p->perturb_coeffs   = this->perturb_coeffs;
    p->scenario_names   = this->scenario_names;
    p->probability      = this->probability;
    p->invest_coeff     = this->invest_coeff;
    p->first_stage_map  = this->first_stage_map;
    p->second_stage_map = this->second_stage_map;
    p->clearDAG();
    if (this->full_model_built)
        p->buildFullModelDAG();
    else
        p->buildDAG();
    return p;
}
