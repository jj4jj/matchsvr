#pragma  once
#include <vector>
#include <cstdint>

class MatchedResult_ST;
class MatchingTeam_ST;
struct MatchingPoolImpl;
struct MatchingPool {
    int         init(int team_member_num);
    int         update(int past_ms, int checked_num = 1000);
    uint32_t    time() const;
    int                         join(const std::vector<uint64_t> & members);
    const MatchingTeam_ST *     quit(uint64_t member_id);
    void                        free_team(const MatchingTeam_ST * team);
    const MatchingTeam_ST *     get_team(size_t team_id);

    const MatchedResult_ST * fetch_results(int * result_num);
    void                     clear_fetched_results();
public:
    static  int   update_matching_object(uint64_t id, int elo, int lv);

private:
    MatchingPoolImpl * impl_;
public:
    MatchingPool();
    ~MatchingPool();
};




