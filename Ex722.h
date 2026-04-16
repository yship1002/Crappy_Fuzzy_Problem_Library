#include "STModel.h"
class Ex722Model:public STModel {
    public:


        Ex722Model(BranchingStrategy branching_strategy);
        Ex722Model(const Ex722Model& other)=default;
        Ex722Model()=default; // default constructor
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