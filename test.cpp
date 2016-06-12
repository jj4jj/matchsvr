

#include "matching.cex.hpp"
#include "matching.h"


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
            mq.update_matching_object(role_id, elo, rand() % 25 + 1);
            members.push_back(role_id);
        }

        if (rand() % 100 < 5){
            //test exit
            //mq.quit(ullNextEntID - 1);
        }

        int ret = mq.join(members);
        if (ret){
            printf("ret join :%d \n", ret);
        }
        mq.update(rand() % 100, rand() % 100);
        if (i % 5 == 0){
            int n = rand() % 20;
            const MatchedResult_ST * r = mq.fetch_results(&n);
            for (int j = 0; j < n; ++j){
                const MatchingTeam_ST * lt = mq.get_team(r[j].teams[0]);
                const MatchingTeam_ST * rt = mq.get_team(r[j].teams[1]);

                printf("Got an pair ! [%d\t%d\tin\t%d]  <->  [%d\t%d\tin\t%d] \n",
                    lt->point.elo, lt->point.level, mq.time() - lt->join_time,
                    rt->point.elo, rt->point.level, mq.time() - rt->join_time);
            }
            mq.clear_fetched_results();
        }
    }
    return 0;
}
