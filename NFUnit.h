#include "STModel.h"

class NFUnit : public STModel {
    public:
        NFUnit(BranchingStrategy branchingStrategy);
        NFUnit(const NFUnit& other) = default;
        NFUnit() = default;

        double probability;
        Ipopt::SmartPtr<STModel> clone() override;
        void buildDAG() override;
        void buildFullModelDAG() override;

        template<class Archive>
        void serialize(Archive& ar) {
            ar(CEREAL_NVP(first_stage_IX), CEREAL_NVP(second_stage_IX));
        }
};