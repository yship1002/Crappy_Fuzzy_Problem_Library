#include "STModel.h"
class TTT:public STModel {
    public:


        TTT(BranchingStrategy branching_strategy);
        TTT(const TTT& other)=default;
        TTT()=default; // default constructor
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