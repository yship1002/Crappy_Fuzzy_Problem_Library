#include "STModel.h"
class Ex844Model:public STModel {
    public:


        Ex844Model(BranchingStrategy branching_strategy);
        Ex844Model(const Ex844Model& other)=default;
        Ex844Model()=default; // default constructor

        std::map<ScenarioNames, double> perturb_a;
        std::map<ScenarioNames, double> perturb_b;
        std::map<ScenarioNames, double> perturb_c;
        std::map<ScenarioNames, double> perturb_d;
        std::map<ScenarioNames, double> perturb_e;
        std::map<ScenarioNames, double> target_shift_x1;
        std::map<ScenarioNames, double> target_shift_x2;
        std::map<ScenarioNames, double> target_shift_x3;
        std::map<ScenarioNames, double> target_shift_x4;
        std::map<ScenarioNames, double> target_shift_x5;
        std::map<ScenarioNames, double> prob;

        /// A vector of ScenarioNames
        double probability;
        Ipopt::SmartPtr<STModel> clone() override;
        void buildDAG() override;
        void buildFullModelDAG() override;

        /**
        * @brief tell cereal what data to save
        * @tparam Archive 
        * @param ar 
        */
        template<class Archive>
        void serialize(Archive& ar) {
        CEREAL_NVP(first_stage_IX), CEREAL_NVP(second_stage_IX);
        }

};