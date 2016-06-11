#pragma  once
#include <vector>
#include <cstdint>



class MatchedResult_ST;
class MatchingTeam_ST;
struct MatchingQueueImpl;
struct MatchingQueue {
    int     init(int team_member_num);
    int     join(const std::vector<uint64_t> & members);
    MatchingTeam_ST *     exit(uint64_t member_id);
    int     update(int past_ms, int checked_num = 1000);
    const MatchedResult_ST * fetch_results(int * result_num);
public:
    static  int   register_matching_object(uint64_t id, int elo, int lv);

private:
    MatchingQueueImpl * impl_;
    MatchingQueue();
    ~MatchingQueue();

};




