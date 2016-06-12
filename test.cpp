

#include "matching.cex.hpp"
#include "matching.h"

static inline int  _elo_to_level(int elo){
    if (elo <= 500){
        return 0;
    }
    if (elo >= 3000){
        return MACTCHING_QUEUE_MAX_LEVEL;
    }
    //todo
    return (elo - 600) / 100;
}
int main(){

    MatchingPool mq;
    mq.init(3);

    static uint64_t ullNextEntID = 10000;
    std::vector<uint64_t>   members;
    for (int i = 0; i < 10000; ++i){
        int iPoint = rand() % 2500 + 500;

        int n = rand() % 3 + 1;
        members.clear();
        for (int j = 0; j < n; ++j){
            uint64_t role_id = ullNextEntID++;
            uint32_t elo = iPoint;
            mq.update_matching_object(role_id, elo, _elo_to_level(elo));// rand() % 25 + 1);
            members.push_back(role_id);
        }
        int ret = mq.join(members);
        if (ret){
            printf("ret join :%d \n", ret);
        }

        if (rand() % 100 < 5){
            //test exit
            const MatchingTeam_ST * team = mq.quit(ullNextEntID - 1);
            assert(team);
            printf("quit team success role:%lu!\n", ullNextEntID - 1);
            mq.free_team(team);
        }

        mq.update(rand() % 100, rand() % 100);
        if (i % 5 == 0){
            size_t n = rand() % 20;            
            const MatchedResult_ST * r = mq.fetch_results(&n);
            printf("fetched result num:%zu\n", n);
            for (size_t j = 0; j < n; ++j){
                printf("process result pair ! %zu <-> %zu\n", 
                    r[j].teams[0],
                    r[j].teams[1]);
                const MatchingTeam_ST * lt = mq.get_team(r[j].teams[0]);
                const MatchingTeam_ST * rt = mq.get_team(r[j].teams[1]);
                printf("process result pair ! %zu#[%d\t%d\tin\t%d]  <->  %zu#[%d\t%d\tin\t%d] \n",
                    r[j].teams[0],
                    lt->point.elo, lt->point.level, mq.time() - lt->join_time,
                    r[j].teams[1],
                    rt->point.elo, rt->point.level, mq.time() - rt->join_time);
            }
            mq.clear_fetched_results();
        }
        if (rand() % 100 < 10){
            //test exit
            const MatchingTeam_ST * team = mq.quit(ullNextEntID - 1);
            if (team){
                printf("quit team success role:%lu!\n", ullNextEntID - 1);
                mq.free_team(team);
            }
            else {
                printf("quit team fail role:%lu!\n", ullNextEntID - 1);
            }
        }


    }
    return 0;
}

