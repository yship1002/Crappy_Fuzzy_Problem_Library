// #include "Ex844.h"

// // ============================================================================
// // Variable index convention (same as template):
// //   X[0]=x6   X[1]=x7   X[2]=x8   X[3]=x9   X[4]=x10  X[5]=x11  X[6]=x12
// //   X[7]=x13  X[8]=x14  X[9]=x15  X[10]=x16 X[11]=x17
// //   X[12]=x1[s] X[13]=x2[s] X[14]=x3[s] X[15]=x4[s] X[16]=x5[s]
// //
// // Constraint mapping from the Pyomo model:
// //   First-stage only (unchanged vs. the old model):
// //     c1 <-> e4   c2 <-> e5   c3 <-> e6
// //     c4 <-> e9   c5 <-> e10  c6 <-> e11  c7 <-> e12
// //   Scenario-indexed (perturb_a..e are now five INDEPENDENT draws, not one
// //   shared scalar; e8 is a new active constraint that didn't exist before):
// //     c8  <-> e1 (perturb_a, hits x1[s])
// //     c9  <-> e2 (perturb_b, hits x3[s])
// //     c10 <-> e3 (perturb_c, hits x5[s])
// //     c11 <-> e7 (perturb_d, hits x2[s])
// //     c12 <-> e8 (perturb_e, hits x4[s])  -- NEW, was commented out before
// // ============================================================================

// Ex844Model::Ex844Model(BranchingStrategy branching_strategy):STModel() {

//     this->branching_strategy = branching_strategy;

//     this->scenario_names = {
//         ScenarioNames::SCENARIO1,  ScenarioNames::SCENARIO2,  ScenarioNames::SCENARIO3,
//         ScenarioNames::SCENARIO4,  ScenarioNames::SCENARIO5,  ScenarioNames::SCENARIO6,
//         ScenarioNames::SCENARIO7,  ScenarioNames::SCENARIO8,  ScenarioNames::SCENARIO9,
//         ScenarioNames::SCENARIO10
//     };
//     this->scenario_name = ScenarioNames::SCENARIO1; // default

//     // ------------------------------------------------------------------
//     // CHANGE 1: five independent per-constraint perturbation vectors,
//     // replacing the single shared "perturb" scalar the old model used
//     // for every scenario constraint.
//     //   perturb_a -> e1 (x1[s])   perturb_b -> e2 (x3[s])
//     //   perturb_c -> e3 (x5[s])   perturb_d -> e7 (x2[s])
//     //   perturb_e -> e8 (x4[s])
//     // ------------------------------------------------------------------
//     this->perturb_a = {
//         {ScenarioNames::SCENARIO1,  0.16435402478574515},
//         {ScenarioNames::SCENARIO2, -0.27625594348335564},
//         {ScenarioNames::SCENARIO3, -0.5508317712765664},
//         {ScenarioNames::SCENARIO4, -0.580166837365765},
//         {ScenarioNames::SCENARIO5,  0.3759242870403269},
//         {ScenarioNames::SCENARIO6,  0.49530669273326605},
//         {ScenarioNames::SCENARIO7,  0.12796293092061584},
//         {ScenarioNames::SCENARIO8,  0.27539587318079806},
//         {ScenarioNames::SCENARIO9,  0.052349989758507444},
//         {ScenarioNames::SCENARIO10, 0.5220869085453218}
//     };
//     this->perturb_b = {
//         {ScenarioNames::SCENARIO1,  0.3790242649458386},
//         {ScenarioNames::SCENARIO2, -0.5967137997958223},
//         {ScenarioNames::SCENARIO3,  0.4288851319050832},
//         {ScenarioNames::SCENARIO4, -0.5596973096334428},
//         {ScenarioNames::SCENARIO5,  0.2755865357159329},
//         {ScenarioNames::SCENARIO6, -0.3892132552769292},
//         {ScenarioNames::SCENARIO7,  0.4358147068198639},
//         {ScenarioNames::SCENARIO8,  0.049753464298910054},
//         {ScenarioNames::SCENARIO9, -0.24034573135513826},
//         {ScenarioNames::SCENARIO10,-0.09277533456280986}
//     };
//     this->perturb_c = {
//         {ScenarioNames::SCENARIO1, -0.5660163946254444},
//         {ScenarioNames::SCENARIO2, -0.45086006820052327},
//         {ScenarioNames::SCENARIO3,  0.2047492976323564},
//         {ScenarioNames::SCENARIO4,  0.1766274138891001},
//         {ScenarioNames::SCENARIO5,  0.13846213377750463},
//         {ScenarioNames::SCENARIO6, -0.13958693488573987},
//         {ScenarioNames::SCENARIO7,  0.5966519229470532},
//         {ScenarioNames::SCENARIO8,  0.577002406531476},
//         {ScenarioNames::SCENARIO9,  0.22265038137683363},
//         {ScenarioNames::SCENARIO10, 0.18055113152137955}
//     };
//     this->perturb_d = {
//         {ScenarioNames::SCENARIO1,  0.22613607668512808},
//         {ScenarioNames::SCENARIO2, -0.13329429122507544},
//         {ScenarioNames::SCENARIO3, -0.4378841939731065},
//         {ScenarioNames::SCENARIO4,  0.26578600823289805},
//         {ScenarioNames::SCENARIO5,  0.030425186970871064},
//         {ScenarioNames::SCENARIO6, -0.2277097493292532},
//         {ScenarioNames::SCENARIO7, -0.016997569401853107},
//         {ScenarioNames::SCENARIO8,  0.4673854012188003},
//         {ScenarioNames::SCENARIO9,  0.5208522191474996},
//         {ScenarioNames::SCENARIO10,-0.1706457639491157}
//     };
//     this->perturb_e = {
//         {ScenarioNames::SCENARIO1,  0.08583579687571312},
//         {ScenarioNames::SCENARIO2, -0.2137567307088694},
//         {ScenarioNames::SCENARIO3,  0.11316003623963615},
//         {ScenarioNames::SCENARIO4, -0.1945065293914401},
//         {ScenarioNames::SCENARIO5, -0.13005719936620652},
//         {ScenarioNames::SCENARIO6,  0.46832922240575076},
//         {ScenarioNames::SCENARIO7, -0.3274108877599443},
//         {ScenarioNames::SCENARIO8,  0.14782457362325088},
//         {ScenarioNames::SCENARIO9, -0.4991815877011382},
//         {ScenarioNames::SCENARIO10, 0.39917297718407735}
//     };

//     // ------------------------------------------------------------------
//     // CHANGE 2: scenario-dependent objective targets. Stored here as
//     // shifts relative to the fixed base target point; the base constants
//     // (base_x1..base_x5) are applied where the objective is built.
//     // ------------------------------------------------------------------
//     this->target_shift_x1 = {
//         {ScenarioNames::SCENARIO1,  0.2739233746429086},
//         {ScenarioNames::SCENARIO2, -0.4604265724722594},
//         {ScenarioNames::SCENARIO3, -0.9180529521276106},
//         {ScenarioNames::SCENARIO4, -0.9669447289429418},
//         {ScenarioNames::SCENARIO5,  0.6265404784005448},
//         {ScenarioNames::SCENARIO6,  0.8255111545554434},
//         {ScenarioNames::SCENARIO7,  0.21327155153435973},
//         {ScenarioNames::SCENARIO8,  0.4589931219679968},
//         {ScenarioNames::SCENARIO9,  0.08724998293084574},
//         {ScenarioNames::SCENARIO10, 0.8701448475755365}
//     };
//     this->target_shift_x2 = {
//         {ScenarioNames::SCENARIO1,  0.6317071082430643},
//         {ScenarioNames::SCENARIO2, -0.9945229996597038},
//         {ScenarioNames::SCENARIO3,  0.7148085531751387},
//         {ScenarioNames::SCENARIO4, -0.9328288493890713},
//         {ScenarioNames::SCENARIO5,  0.45931089285988813},
//         {ScenarioNames::SCENARIO6, -0.648688758794882},
//         {ScenarioNames::SCENARIO7,  0.7263578446997732},
//         {ScenarioNames::SCENARIO8,  0.08292244049818343},
//         {ScenarioNames::SCENARIO9, -0.40057621892523043},
//         {ScenarioNames::SCENARIO10,-0.1546255576046831}
//     };
//     this->target_shift_x3 = {
//         {ScenarioNames::SCENARIO1, -0.9433606577090741},
//         {ScenarioNames::SCENARIO2, -0.7514334470008721},
//         {ScenarioNames::SCENARIO3,  0.34124882938726064},
//         {ScenarioNames::SCENARIO4,  0.2943790231485002},
//         {ScenarioNames::SCENARIO5,  0.2307702229625077},
//         {ScenarioNames::SCENARIO6, -0.23264489147623313},
//         {ScenarioNames::SCENARIO7,  0.994419871578422},
//         {ScenarioNames::SCENARIO8,  0.9616706775524602},
//         {ScenarioNames::SCENARIO9,  0.3710839689613894},
//         {ScenarioNames::SCENARIO10, 0.3009185525356326}
//     };
//     this->target_shift_x4 = {
//         {ScenarioNames::SCENARIO1,  0.37689346114188016},
//         {ScenarioNames::SCENARIO2, -0.22215715204179243},
//         {ScenarioNames::SCENARIO3, -0.7298069899551776},
//         {ScenarioNames::SCENARIO4,  0.4429766803881634},
//         {ScenarioNames::SCENARIO5,  0.050708644951451776},
//         {ScenarioNames::SCENARIO6, -0.3795162488820887},
//         {ScenarioNames::SCENARIO7, -0.028329282336421846},
//         {ScenarioNames::SCENARIO8,  0.7789756686980005},
//         {ScenarioNames::SCENARIO9,  0.8680870319124994},
//         {ScenarioNames::SCENARIO10,-0.28440960658185954}
//     };
//     this->target_shift_x5 = {
//         {ScenarioNames::SCENARIO1,  0.14305966145952187},
//         {ScenarioNames::SCENARIO2, -0.3562612178481157},
//         {ScenarioNames::SCENARIO3,  0.1886000603993936},
//         {ScenarioNames::SCENARIO4, -0.3241775489857335},
//         {ScenarioNames::SCENARIO5, -0.21676199894367754},
//         {ScenarioNames::SCENARIO6,  0.7805487040095846},
//         {ScenarioNames::SCENARIO7, -0.5456848129332406},
//         {ScenarioNames::SCENARIO8,  0.24637428937208483},
//         {ScenarioNames::SCENARIO9, -0.8319693128352303},
//         {ScenarioNames::SCENARIO10, 0.6652882953067956}
//     };

//     // ------------------------------------------------------------------
//     // CHANGE 3: non-uniform scenario probabilities, replacing the old
//     // scalar this->probability (which assumed equal weight per scenario).
//     // Still sums to 1.
//     // ------------------------------------------------------------------
//     this->prob = {
//         {ScenarioNames::SCENARIO1,  0.14337344},
//         {ScenarioNames::SCENARIO2,  0.08236041},
//         {ScenarioNames::SCENARIO3,  0.15333038},
//         {ScenarioNames::SCENARIO4,  0.06222044},
//         {ScenarioNames::SCENARIO5,  0.09313739},
//         {ScenarioNames::SCENARIO6,  0.07243643},
//         {ScenarioNames::SCENARIO7,  0.10586093},
//         {ScenarioNames::SCENARIO8,  0.14440115},
//         {ScenarioNames::SCENARIO9,  0.08138826},
//         {ScenarioNames::SCENARIO10, 0.06149118}
//     };

//     // First-stage variable bounds (unchanged)
//     this->first_stage_IX = {
//         mc::Interval(-2, 0),
//         mc::Interval(0.5, 2.5),
//         mc::Interval(-1.5, 0.5),
//         mc::Interval(0.2, 2.2),
//         mc::Interval(-1.2, 0.8),
//         mc::Interval(0.1, 2.1),
//         mc::Interval(-1.1, 0.9),
//         mc::Interval(0, 1),
//         mc::Interval(0, 1),
//         mc::Interval(1.1, 1.3),
//         mc::Interval(0, 1),
//         mc::Interval(0, 1)
//     };

//     // Second-stage (recourse) variable bounds, per scenario (unchanged)
//     this->second_stage_IX = {
//         mc::Interval(4, 6),
//         mc::Interval(-6, -4),
//         mc::Interval(2, 4),
//         mc::Interval(-3, -1),
//         mc::Interval(1, 3),
//     };
// };

// void Ex844Model::buildDAG() {

//     // Fixed base target point; scenario-specific offsets are added below.
//     const double base_x1 = -5.0, base_x2 = 5.0, base_x3 = -3.0, base_x4 = 2.0, base_x5 = -2.0;

//     for (auto& scenario_name : this->scenario_names) {
//         int n_first_stage_vars = this->first_stage_IX.size();
//         const int nvars = n_first_stage_vars + this->second_stage_IX.size();

//         this->X[scenario_name].resize(nvars);
//         for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);

//         // NOTE: fixed vs. the template, which incorrectly read
//         // this->perturb[this->scenario_name] (always the default scenario)
//         // instead of the current loop scenario. Each scenario now pulls
//         // its own five independent perturbation values.
//         double p_a = this->perturb_a[scenario_name];
//         double p_b = this->perturb_b[scenario_name];
//         double p_c = this->perturb_c[scenario_name];
//         double p_d = this->perturb_d[scenario_name];
//         double p_e = this->perturb_e[scenario_name];

//         double p_s = this->prob[scenario_name];
//         double ts1 = this->target_shift_x1[scenario_name];
//         double ts2 = this->target_shift_x2[scenario_name];
//         double ts3 = this->target_shift_x3[scenario_name];
//         double ts4 = this->target_shift_x4[scenario_name];
//         double ts5 = this->target_shift_x5[scenario_name];

//         // Indices: X[0]=x6 ... X[11]=x17, X[12]=x1[s] ... X[16]=x5[s]
//         mc::FFVar c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;
//         mc::FFVar nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10, nc11, nc12;

//         double scale = 1;

//         // -------- first-stage-only constraints (unchanged: e4,e5,e6,e9,e10,e11,e12) --------
//         c1 = this->X[scenario_name][8] - pow(0.628318, this->X[scenario_name][9]) * (this->X[scenario_name][1] - this->X[scenario_name][7]);
//         c1 = scale * c1;
//         nc1 = -c1;

//         c2 = this->X[scenario_name][8] - pow(0.7853975, this->X[scenario_name][9]) * (this->X[scenario_name][3] - this->X[scenario_name][7]);
//         c2 = scale * c2;
//         nc2 = -c2;

//         c3 = this->X[scenario_name][8] - pow(0.942477, this->X[scenario_name][9]) * (this->X[scenario_name][5] - this->X[scenario_name][7]);
//         c3 = scale * c3;
//         nc3 = -c3;

//         c4 = -this->X[scenario_name][11] - pow(0.4712385, this->X[scenario_name][9]) * (this->X[scenario_name][0] - 0.4712385 * this->X[scenario_name][10]);
//         c4 = scale * c4;
//         nc4 = -c4;

//         c5 = -this->X[scenario_name][11] - pow(0.628318, this->X[scenario_name][9]) * (this->X[scenario_name][2] - 0.628318 * this->X[scenario_name][10]);
//         c5 = scale * c5;
//         nc5 = -c5;

//         c6 = -this->X[scenario_name][11] - pow(0.7853975, this->X[scenario_name][9]) * (this->X[scenario_name][4] - 0.7853975 * this->X[scenario_name][10]);
//         c6 = scale * c6;
//         nc6 = -c6;

//         c7 = -this->X[scenario_name][11] - pow(0.942477, this->X[scenario_name][9]) * (this->X[scenario_name][6] - 0.942477 * this->X[scenario_name][10]);
//         c7 = scale * c7;
//         nc7 = -c7;

//         // -------- scenario constraints (e1,e2,e3,e7,e8), each with its own perturbation --------
//         c8 = this->X[scenario_name][8] - pow(0.1570795, this->X[scenario_name][9]) * (this->X[scenario_name][12] - this->X[scenario_name][7] + p_a);
//         c8 = scale * c8;
//         nc8 = -c8;

//         c9 = this->X[scenario_name][8] - pow(0.314159, this->X[scenario_name][9]) * (this->X[scenario_name][14] - this->X[scenario_name][7] + p_b);
//         c9 = scale * c9;
//         nc9 = -c9;

//         c10 = this->X[scenario_name][8] - pow(0.4712385, this->X[scenario_name][9]) * (this->X[scenario_name][16] - this->X[scenario_name][7] + p_c);
//         c10 = scale * c10;
//         nc10 = -c10;

//         c11 = -this->X[scenario_name][11] - pow(0.1570795, this->X[scenario_name][9]) * (this->X[scenario_name][13] - 0.1570795 * this->X[scenario_name][10] + p_d);
//         c11 = scale * c11;
//         nc11 = -c11;

//         // NEW: e8 was commented out (x4[s] unconstrained) in the old model; now active.
//         c12 = -this->X[scenario_name][11] - pow(0.314159, this->X[scenario_name][9]) * (this->X[scenario_name][15] - 0.314159 * this->X[scenario_name][10] + p_e);
//         c12 = scale * c12;
//         nc12 = -c12;

//         // -------- objective: scenario-specific targets and scenario-specific probability --------
//         mc::FFVar objective = 1000*p_s * (
//               pow((base_x1 + ts1) + this->X[scenario_name][12], 2)
//             + pow((base_x2 + ts2) + this->X[scenario_name][13], 2)
//             + pow((base_x3 + ts3) + this->X[scenario_name][14], 2)
//             + pow((base_x4 + ts4) + this->X[scenario_name][15], 2)
//             + pow((base_x5 + ts5) + this->X[scenario_name][16], 2)
//             + pow(1 + this->X[scenario_name][0], 2)
//             + pow(-1.5 + this->X[scenario_name][1], 2)
//             + pow(0.5 + this->X[scenario_name][2], 2)
//             + pow(-1.2 + this->X[scenario_name][3], 2)
//             + pow(0.2 + this->X[scenario_name][4], 2)
//             + pow(-1.1 + this->X[scenario_name][5], 2)
//             + pow(0.1 + this->X[scenario_name][6], 2)
//         );

//         this->F[scenario_name] = {
//             objective,
//             c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
//             nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10, nc11, nc12
//         };
//     }
// }

// void Ex844Model::buildFullModelDAG() {

//     const double base_x1 = -5.0, base_x2 = 5.0, base_x3 = -3.0, base_x4 = 2.0, base_x5 = -2.0;

//     int n_first_stage_vars = this->first_stage_IX.size();
//     // NOTE: fixed vs. the template. second_stage_IX only ever holds bounds
//     // for ONE scenario's 5 recourse vars (reused per scenario), so the
//     // per-scenario count is just its size -- dividing by scenario_names.size()
//     // (as the template did) truncates to 0 under integer division and
//     // breaks the model.
//     int n_second_stage_vars = this->second_stage_IX.size()/this->scenario_names.size();
//     int n_scenarios = this->scenario_names.size();
//     int nvars = n_first_stage_vars +   n_scenarios*n_second_stage_vars;

//     this->X[ScenarioNames::SCENARIO1].resize(nvars);

    
//     for (int i = 0; i < nvars; ++i) this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);

//     mc::FFVar objective = 0;

//     for (int s_idx = 0; s_idx < n_scenarios; ++s_idx) {
//         ScenarioNames scenario_name = this->scenario_names[s_idx];
//         int second_stage_start_idx = n_first_stage_vars + s_idx * n_second_stage_vars;

//         double p_a = this->perturb_a[scenario_name];
//         double p_b = this->perturb_b[scenario_name];
//         double p_c = this->perturb_c[scenario_name];
//         double p_d = this->perturb_d[scenario_name];
//         double p_e = this->perturb_e[scenario_name];

//         double p_s = this->prob[scenario_name];
//         double ts1 = this->target_shift_x1[scenario_name];
//         double ts2 = this->target_shift_x2[scenario_name];
//         double ts3 = this->target_shift_x3[scenario_name];
//         double ts4 = this->target_shift_x4[scenario_name];
//         double ts5 = this->target_shift_x5[scenario_name];

//         auto& X = this->X[ScenarioNames::SCENARIO1];

//         mc::FFVar c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;
//         mc::FFVar nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10, nc11, nc12;

//         double scale = 1;

//         c1 = X[8] - pow(0.628318, X[9]) * (X[1] - X[7]);
//         c1 = scale * c1;
//         nc1 = -c1;

//         c2 = X[8] - pow(0.7853975, X[9]) * (X[3] - X[7]);
//         c2 = scale * c2;
//         nc2 = -c2;

//         c3 = X[8] - pow(0.942477, X[9]) * (X[5] - X[7]);
//         c3 = scale * c3;
//         nc3 = -c3;

//         c4 = -X[11] - pow(0.4712385, X[9]) * (X[0] - 0.4712385 * X[10]);
//         c4 = scale * c4;
//         nc4 = -c4;

//         c5 = -X[11] - pow(0.628318, X[9]) * (X[2] - 0.628318 * X[10]);
//         c5 = scale * c5;
//         nc5 = -c5;

//         c6 = -X[11] - pow(0.7853975, X[9]) * (X[4] - 0.7853975 * X[10]);
//         c6 = scale * c6;
//         nc6 = -c6;

//         c7 = -X[11] - pow(0.942477, X[9]) * (X[6] - 0.942477 * X[10]);
//         c7 = scale * c7;
//         nc7 = -c7;

//         c8 = X[8] - pow(0.1570795, X[9]) * (X[second_stage_start_idx] - X[7] + p_a);
//         c8 = scale * c8;
//         nc8 = -c8;

//         c9 = X[8] - pow(0.314159, X[9]) * (X[second_stage_start_idx + 2] - X[7] + p_b);
//         c9 = scale * c9;
//         nc9 = -c9;

//         c10 = X[8] - pow(0.4712385, X[9]) * (X[second_stage_start_idx + 4] - X[7] + p_c);
//         c10 = scale * c10;
//         nc10 = -c10;

//         c11 = -X[11] - pow(0.1570795, X[9]) * (X[second_stage_start_idx + 1] - 0.1570795 * X[10] + p_d);
//         c11 = scale * c11;
//         nc11 = -c11;

//         // NEW: e8, previously inactive
//         c12 = -X[11] - pow(0.314159, X[9]) * (X[second_stage_start_idx + 3] - 0.314159 * X[10] + p_e);
//         c12 = scale * c12;
//         nc12 = -c12;

//         objective += 1000*p_s * (
//               pow((base_x1 + ts1) + X[second_stage_start_idx], 2)
//             + pow((base_x2 + ts2) + X[second_stage_start_idx + 1], 2)
//             + pow((base_x3 + ts3) + X[second_stage_start_idx + 2], 2)
//             + pow((base_x4 + ts4) + X[second_stage_start_idx + 3], 2)
//             + pow((base_x5 + ts5) + X[second_stage_start_idx + 4], 2)
//             + pow(1 + X[0], 2)
//             + pow(-1.5 + X[1], 2)
//             + pow(0.5 + X[2], 2)
//             + pow(-1.2 + X[3], 2)
//             + pow(0.2 + X[4], 2)
//             + pow(-1.1 + X[5], 2)
//             + pow(0.1 + X[6], 2)
//         );

//         std::vector<mc::FFVar> scenario_constraints = {
//             c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
//             nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10, nc11, nc12
//         };
//         this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].end(), scenario_constraints.begin(), scenario_constraints.end());
//     }

//     this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].begin(), objective);
//     this->full_model_built = true;
// }

// Ipopt::SmartPtr<STModel> Ex844Model::clone() {
//     Ipopt::SmartPtr<Ex844Model> p = new Ex844Model();

//     p->scenario_name = this->scenario_name;
//     p->first_stage_IX = this->first_stage_IX;
//     p->second_stage_IX = this->second_stage_IX;
//     p->scenario_names = this->scenario_names;

//     p->perturb_a = this->perturb_a;
//     p->perturb_b = this->perturb_b;
//     p->perturb_c = this->perturb_c;
//     p->perturb_d = this->perturb_d;
//     p->perturb_e = this->perturb_e;

//     p->prob = this->prob;
//     p->target_shift_x1 = this->target_shift_x1;
//     p->target_shift_x2 = this->target_shift_x2;
//     p->target_shift_x3 = this->target_shift_x3;
//     p->target_shift_x4 = this->target_shift_x4;
//     p->target_shift_x5 = this->target_shift_x5;

//     p->clearDAG(); // clear the DAG of the cloned model
//     if (this->full_model_built) {
//         p->buildFullModelDAG();
//     } else {
//         p->buildDAG();
//     }
//     return p;
// }

#include "Ex844.h"

// ============================================================================
// SIMPLIFIED two-stage nonconvex test problem (for diagnosing decomposition)
//
// CONVENTION ASSUMED: every entry of F after the objective is a constraint
//                     g(x) <= 0.  (If your framework uses g(x) >= 0, flip signs.)
//
// Variables (per scenario subproblem):
//   First stage (3):   X[0]=z1 in [0,2]   X[1]=z2 in [0,2]   X[2]=z3 in [1,2]
//   Second stage (2):  X[3]=x1[s] in [0,4]   X[4]=x2[s] in [0,4]
//
// Model:
//   min  sum_s 1000*p_s * [ (x1[s]-T1[s])^2 + (x2[s]-T2[s])^2
//                           + (z1-1)^2 + (z2-1)^2 + (z3-1.5)^2 ]
//   s.t. g0:  0.5 - z1*z2                                   <= 0   (first stage, bilinear)
//        g1:  z1*x1[s] + x2[s] - (3 + pa[s])               <= 0   (bilinear coupling)
//        g2:  z2*x2[s] - 0.5^z3 * x1[s] - (1 + pb[s])      <= 0   (bilinear + exponential)
//   with T1[s] = 2.5 + ts1[s],  T2[s] = 2.0 + ts2[s].
//
// Why it is always feasible:
//   - g0 holds at z=(1,1,*).
//   - For ANY first-stage point, x1[s]=x2[s]=0 satisfies g1 and g2 strictly in
//     every scenario (|pa|,|pb| < 0.6), so the recourse is relatively complete.
//   - Known strictly feasible point for the full model: z=(1,1,1.5), x[s]=(0,0).
//
// Why the recourse is a real optimization:
//   For fixed z, each scenario is a 2-variable QP with 2 active-able linear
//   constraints -- x1[s], x2[s] are chosen, not determined. Nonconvexity lives
//   in the first/second-stage coupling (z*x bilinear terms, 0.5^z3 * x1).
//
// Reference solution (see verify_ex844.py; exact recourse on an 81^3 first-stage
// grid + local refinement, cross-checked with 150-start NLP on the extensive form):
//   GLOBAL full-model optimum  = 583.1556  at z* = (0.4783, 1.0453, 1.4876)
//   Sum of independent scenario optima (wait-and-see bound) = 500.7171
//   Gap = 14.1%
// ============================================================================

namespace {
    const double T1_BASE = 2.5, T2_BASE = 2.0;
}

Ex844Model::Ex844Model(BranchingStrategy branching_strategy):STModel() {

    this->branching_strategy = branching_strategy;

    this->scenario_names = {
        ScenarioNames::SCENARIO1,  ScenarioNames::SCENARIO2,  ScenarioNames::SCENARIO3,
        ScenarioNames::SCENARIO4,  ScenarioNames::SCENARIO5,  ScenarioNames::SCENARIO6,
        ScenarioNames::SCENARIO7,  ScenarioNames::SCENARIO8,  ScenarioNames::SCENARIO9,
        ScenarioNames::SCENARIO10
    };
    this->scenario_name = ScenarioNames::SCENARIO1;

    // Constraint RHS perturbations (same numbers as before, only two are used now)
    this->perturb_a = {
        {ScenarioNames::SCENARIO1,  0.16435402478574515},
        {ScenarioNames::SCENARIO2, -0.27625594348335564},
        {ScenarioNames::SCENARIO3, -0.5508317712765664},
        {ScenarioNames::SCENARIO4, -0.580166837365765},
        {ScenarioNames::SCENARIO5,  0.3759242870403269},
        {ScenarioNames::SCENARIO6,  0.49530669273326605},
        {ScenarioNames::SCENARIO7,  0.12796293092061584},
        {ScenarioNames::SCENARIO8,  0.27539587318079806},
        {ScenarioNames::SCENARIO9,  0.052349989758507444},
        {ScenarioNames::SCENARIO10, 0.5220869085453218}
    };
    this->perturb_b = {
        {ScenarioNames::SCENARIO1,  0.3790242649458386},
        {ScenarioNames::SCENARIO2, -0.5967137997958223},
        {ScenarioNames::SCENARIO3,  0.4288851319050832},
        {ScenarioNames::SCENARIO4, -0.5596973096334428},
        {ScenarioNames::SCENARIO5,  0.2755865357159329},
        {ScenarioNames::SCENARIO6, -0.3892132552769292},
        {ScenarioNames::SCENARIO7,  0.4358147068198639},
        {ScenarioNames::SCENARIO8,  0.049753464298910054},
        {ScenarioNames::SCENARIO9, -0.24034573135513826},
        {ScenarioNames::SCENARIO10,-0.09277533456280986}
    };

    // Objective target shifts (same numbers as before, only two are used now)
    this->target_shift_x1 = {
        {ScenarioNames::SCENARIO1,  0.2739233746429086},
        {ScenarioNames::SCENARIO2, -0.4604265724722594},
        {ScenarioNames::SCENARIO3, -0.9180529521276106},
        {ScenarioNames::SCENARIO4, -0.9669447289429418},
        {ScenarioNames::SCENARIO5,  0.6265404784005448},
        {ScenarioNames::SCENARIO6,  0.8255111545554434},
        {ScenarioNames::SCENARIO7,  0.21327155153435973},
        {ScenarioNames::SCENARIO8,  0.4589931219679968},
        {ScenarioNames::SCENARIO9,  0.08724998293084574},
        {ScenarioNames::SCENARIO10, 0.8701448475755365}
    };
    this->target_shift_x2 = {
        {ScenarioNames::SCENARIO1,  0.6317071082430643},
        {ScenarioNames::SCENARIO2, -0.9945229996597038},
        {ScenarioNames::SCENARIO3,  0.7148085531751387},
        {ScenarioNames::SCENARIO4, -0.9328288493890713},
        {ScenarioNames::SCENARIO5,  0.45931089285988813},
        {ScenarioNames::SCENARIO6, -0.648688758794882},
        {ScenarioNames::SCENARIO7,  0.7263578446997732},
        {ScenarioNames::SCENARIO8,  0.08292244049818343},
        {ScenarioNames::SCENARIO9, -0.40057621892523043},
        {ScenarioNames::SCENARIO10,-0.1546255576046831}
    };

    // Scenario probabilities (sum to 1)
    this->prob = {
        {ScenarioNames::SCENARIO1,  0.14337344},
        {ScenarioNames::SCENARIO2,  0.08236041},
        {ScenarioNames::SCENARIO3,  0.15333038},
        {ScenarioNames::SCENARIO4,  0.06222044},
        {ScenarioNames::SCENARIO5,  0.09313739},
        {ScenarioNames::SCENARIO6,  0.07243643},
        {ScenarioNames::SCENARIO7,  0.10586093},
        {ScenarioNames::SCENARIO8,  0.14440115},
        {ScenarioNames::SCENARIO9,  0.08138826},
        {ScenarioNames::SCENARIO10, 0.06149118}
    };

    this->first_stage_IX = {
        mc::Interval(0, 2),   // z1
        mc::Interval(0, 2),   // z2
        mc::Interval(1, 2)    // z3
    };

    // Bounds of ONE scenario's recourse block (reused for every scenario)
    this->second_stage_IX = {
        mc::Interval(0, 4),   // x1[s]
        mc::Interval(0, 4)    // x2[s]
    };
};

void Ex844Model::buildDAG() {

    const int n1 = this->first_stage_IX.size();                 // 3
    const int nvars = n1 + this->second_stage_IX.size();        // 5

    for (auto& scenario_name : this->scenario_names) {

        this->X[scenario_name].resize(nvars);
        for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);
        auto& X = this->X[scenario_name];

        const double pa  = this->perturb_a[scenario_name];
        const double pb  = this->perturb_b[scenario_name];
        const double p_s = this->prob[scenario_name];
        const double T1  = T1_BASE + this->target_shift_x1[scenario_name];
        const double T2  = T2_BASE + this->target_shift_x2[scenario_name];

        const mc::FFVar& z1 = X[0];
        const mc::FFVar& z2 = X[1];
        const mc::FFVar& z3 = X[2];
        const mc::FFVar& x1 = X[n1 + 0];
        const mc::FFVar& x2 = X[n1 + 1];

        mc::FFVar g0 = 0.5 - z1 * z2;
        mc::FFVar g1 = z1 * x1 + x2 - (3.0 + pa);
        mc::FFVar g2 = z2 * x2 - pow(0.5, z3) * x1 - (1.0 + pb);

        mc::FFVar objective = 1000 * p_s * (
              pow(x1 - T1, 2) + pow(x2 - T2, 2)
            + pow(z1 - 1.0, 2) + pow(z2 - 1.0, 2) + pow(z3 - 1.5, 2)
        );

        this->F[scenario_name] = { objective, g0, g1, g2 };
    }
}

void Ex844Model::buildFullModelDAG() {

    const int n1 = this->first_stage_IX.size();                 // 3
    const int nS = this->scenario_names.size();                 // 10
    const int n2 = this->second_stage_IX.size()/nS;                // 2 per scenario
    const int nvars = n1 + nS * n2;                             // 23

    auto& X = this->X[ScenarioNames::SCENARIO1];
    X.resize(nvars);
    for (int i = 0; i < nvars; ++i) X[i].set(&this->DAG[ScenarioNames::SCENARIO1]);

    const mc::FFVar& z1 = X[0];
    const mc::FFVar& z2 = X[1];
    const mc::FFVar& z3 = X[2];

    mc::FFVar objective = 0;
    std::vector<mc::FFVar> G;
    G.push_back(0.5 - z1 * z2);                                 // first-stage constraint, once

    for (int s_idx = 0; s_idx < nS; ++s_idx) {
        ScenarioNames sn = this->scenario_names[s_idx];
        const int b = n1 + s_idx * n2;

        const double pa  = this->perturb_a[sn];
        const double pb  = this->perturb_b[sn];
        const double p_s = this->prob[sn];
        const double T1  = T1_BASE + this->target_shift_x1[sn];
        const double T2  = T2_BASE + this->target_shift_x2[sn];

        const mc::FFVar& x1 = X[b + 0];
        const mc::FFVar& x2 = X[b + 1];

        G.push_back(z1 * x1 + x2 - (3.0 + pa));
        G.push_back(z2 * x2 - pow(0.5, z3) * x1 - (1.0 + pb));

        objective += 1000 * p_s * (
              pow(x1 - T1, 2) + pow(x2 - T2, 2)
            + pow(z1 - 1.0, 2) + pow(z2 - 1.0, 2) + pow(z3 - 1.5, 2)
        );
    }

    this->F[ScenarioNames::SCENARIO1].clear();
    this->F[ScenarioNames::SCENARIO1].push_back(objective);
    this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].end(), G.begin(), G.end());
    this->full_model_built = true;
}

Ipopt::SmartPtr<STModel> Ex844Model::clone() {
    Ipopt::SmartPtr<Ex844Model> p = new Ex844Model(this->branching_strategy);

    p->scenario_name   = this->scenario_name;
    p->scenario_names  = this->scenario_names;
    p->first_stage_IX  = this->first_stage_IX;
    p->second_stage_IX = this->second_stage_IX;

    p->perturb_a       = this->perturb_a;
    p->perturb_b       = this->perturb_b;
    p->target_shift_x1 = this->target_shift_x1;
    p->target_shift_x2 = this->target_shift_x2;
    p->prob            = this->prob;

    p->clearDAG();
    if (this->full_model_built) p->buildFullModelDAG();
    else                        p->buildDAG();
    return p;
}