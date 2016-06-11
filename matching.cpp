#include <time.h>
#include <unordered_map>
#include "matching.cex.hpp"
#include "matching.h"

struct MatchingObjectHash {
    size_t operator ()(const MatchingObject_ST & mo) {
        return mo.id;
    }
};
static struct MatchingEnv {
    pbdcex::mmpool_t<MatchingTeam_ST, MATCHING_QUEUE_MAX_OBJECT_NUM>
        matching_team_pool;
    pbdcex::hashtable_t<MatchingObject_ST, MATCHING_QUEUE_MAX_OBJECT_NUM, MatchingObjectHash>
        matching_object_pool;
    int     inited;
    std::unordered_map<uint64_t, size_t>    member_id_2_team_id;
    MatchingEnv():inited(false){
    }
} s_ENV;

static inline void matching_service_init(){
    if (!s_ENV.inited){
        s_ENV.matching_object_pool.construct();
        s_ENV.matching_team_pool.construct();
    }
}

struct MatchingQueueImpl {
    MatchingQueue_ST queue;
};
MatchingQueue::MatchingQueue(){
    impl_ = new MatchingQueueImpl();
}

MatchingQueue::~MatchingQueue(){
    if (impl_){
        delete impl_;
        impl_ = NULL;
    }
}
int     MatchingQueue::register_matching_object(uint64_t id, int elo, int lv){
    MatchingObject_ST mo;
    mo.construct();
    mo.id = id;
    mo.point.elo = elo;
    mo.point.level = lv;
    if (s_ENV.matching_object_pool.insert(mo)){
        return 0;
    }
    return -1;
}
static inline size_t find_team_id_by_member_id(uint64_t id){

}
#define MQ  (impl_->queue)
int     MatchingQueue::init(int team_member_num){
    matching_service_init();
    //////////////////////////////////////////////
    impl_->queue.construct();
    impl_->queue.team_member_num = team_member_num;
    impl_->queue.cur_time = time(NULL);
    return 0;
}
static inline void _elo_team_init(MatchingTeam_ST & t){
    switch (t.members.count){
    case 1:
        break;
    case 2:
        break;
    case 3:
        break;
    default:
        break;
    }
    //todo get member

}
static inline int  _elo_to_level(int elo){
    if (elo <= 500){
        return 0;
    }
    if (elo >= 3000){
        return MACTCHING_QUEUE_MAX_LEVEL - 1;
    }
    //todo
    return (elo - 600) / 100;
}
static inline bool CheckTeamWaitingCanMerged(const MatchingTeam_ST & t1, const MatchingQueueMatcher_ST & t2){
    int iWaitTime = std::max(MQ.join_ime - t1.join_time,
        MQ.join_ime - t2.join_time);
    int iLvDiff = t1.point.level - t2.m_iLevel;
    if (iLvDiff < 0){ iLvDiff = -iLvDiff; }
    if (iWaitTime >= iLvDiff * 15){
        return true;
    }
    return false;
}
static inline bool _check_result_matched(uint32_t cur_time, MatchingQueueMatcher_ST results[]){
    MatchingTeam_ST * lt = s_ENV.matching_team_pool.ptr(results[MATCHED_RESULT_TEAM_L].team_id);
    MatchingTeam_ST * rt = s_ENV.matching_team_pool.ptr(results[MATCHED_RESULT_TEAM_R].team_id);
    assert(lt && rt);
    int llv = _elo_to_level(results[MATCHED_RESULT_TEAM_L].team_elo);
    int rlv = _elo_to_level(results[MATCHED_RESULT_TEAM_R].team_elo);

    int iWaitTime = std::max(cur_time - lt->join_time, cur_time - rt->join_time);
    int iLvDiff = llv - rlv;
    if (iLvDiff < 0){ iLvDiff = -iLvDiff; }
    if (iWaitTime >= iLvDiff * 15){
        return true;
    }
    return false;
}
int     MatchingQueue::join(const std::vector<uint64_t> & members){
    //ADD TEAM
    if (members.size() > MQ.team_member_num &&
        members.size() == 0){
        return -1;
    }
    MatchingTeam_ST mt;
    mt.construct();
    mt.join_time = MQ.cur_time;
    for (size_t i = 0; i < members.size(); ++i){
        mt.members.lappend(members[i]);
    }
    _elo_team_init(mt);
    size_t mtid = s_ENV.matching_team_pool.alloc();
    if (!mtid){
        return -2;
    }
    *s_ENV.matching_team_pool.ptr(mtid) = mt;

    int iMatchLv = mt.point.level;
    MatchingBucket_ST & bucket = MQ.buckets.list[iMatchLv];

    MatchingQueueMatcher_ST mqst;
    mqst.construct();
    mqst.team_elo = mt.point.elo;
    mqst.team_id = mtid;

    printf("join matching team %d#[%d %d]\n", mt.members.count, iMatchLv, mt.point.elo);
    if (mt.members.count == MQ.team_member_num){
        //insert into matchers
        if (bucket.matchers.count >= MATCHING_TEAM_MAX_NUM_IN_LEVEL){
            s_ENV.matching_team_pool.free(mtid);
            return -3;
        }
        if (bucket.matchers.binsert(mqst)){
            s_ENV.matching_team_pool.free(mtid);
            return -4;
        }
        for (size_t x = 0; x < mt.members.count; ++x){
            s_ENV.member_id_2_team_id[mt.members.list[x]] = mtid;
        }
    }
    else {//waiting merging
        if (bucket.waitings[mt.members.count - 1].teams.count == MACTCHING_BUCKET_MAX_WAITING_TEAM_NUM){
            s_ENV.matching_team_pool.free(mtid);
            return -5; //waiting queue is also full
        }
        if (bucket.waitings[mt.members.count - 1].teams.binsert(mqst)){
            s_ENV.matching_team_pool.free(mtid);
            return -6;
        }
        for (size_t x = 0; x < mt.members.count; ++x){
            s_ENV.member_id_2_team_id[mt.members.list[x]] = mtid;
        }
    }
    //merging
    MatchingTeam_ST merged;
    merged.construct();
    size_t merged_buckets_num = 0;
    struct {
        int level;
        int num_idx;
        int pos_idx;
    } merged_bucket_idxs[MATCHING_TEAM_MAX_NUM_IN_LEVEL];
    //check merging all
    for (int i = MACTCHING_QUEUE_MAX_LEVEL; i >= 0; --i){
        MatchingBucket_ST &    mqbucket = MQ.buckets.list[i];
        for (int j = MQ.team_member_num - 1; j > 0; --j){ //member number max -> min
            int nidx = j - 1;
            MatchingWaitigTeam_ST & wmt = mqbucket.waitings.list[nidx];
            for (int k = wmt.teams.count - 1; k >= 0; --k){
                //team can join  check ?
                if (merged_buckets_num > 0 && !CheckTeamWaitingCanMerged(merged, wmt.teams.list[k])){
                    continue;
                }
                if (merged.members.count + j <= MQ.team_member_num){
                    //merge team
                    merged_bucket_idxs[merged_buckets_num].level = i;
                    merged_bucket_idxs[merged_buckets_num].num_idx = nidx;
                    merged_bucket_idxs[merged_buckets_num].pos_idx = k;
                    ++merged_buckets_num;
                    ///////////////////////////////////////////////////////
                    if (merged.members.count == 0){
                        merged.point.elo = wmt.teams.list[k].team_elo;
                        merged.members.count = j;
                    }
                    else {
                        merged.point.elo = sqrt(merged.point.elo * wmt.teams.list[k].team_elo);
                        merged.members.count += j;
                    }
                }
                //found one
                if (merged.members.count == MQ.team_member_num){
                    printf("found one merge matching team with %d teams [%d]\n",
                        merged_buckets_num, merged.point.elo);
                    int elo_level = _elo_to_level(merged.point.elo);
                    MatchingBucket_ST & merge_bucket = MQ.buckets.list[elo_level];
                    if (merge_bucket.matchers.count >= MATCHING_TEAM_MAX_NUM_IN_LEVEL){
                        memset(&merged, 0, sizeof(merged));
                        merged_buckets_num = 0;
                        continue;
                    }

                    MatchingQueueMatcher_ST nmqst;
                    nmqst.construct();
                    nmqst.team_elo = merged.point.elo;
                    MatchingTeam_ST * nmerging_team = NULL;
                    printf("clear merged \n");
                    //clear merged
                    for (size_t x = 0; x < merged_buckets_num; ++x){
                        MatchingBucket_ST & clear_bucket = MQ.buckets.list[merged_bucket_idxs[x].level];
                        size_t clear_team_id = clear_bucket.waitings.list[merged_bucket_idxs[x].num_idx].\
                            teams.list[merged_bucket_idxs[x].pos_idx].team_id;
                        MatchingTeam_ST *  clear_team = s_ENV.matching_team_pool.ptr(clear_team_id);
                        if (x == 0){
                            nmqst.team_id = clear_team_id;
                            nmerging_team = clear_team;
                            nmerging_team->join_time = MQ.cur_time;
                            nmerging_team->point = merged.point;
                        }
                        else {
                            memcpy(nmerging_team->members.list + nmerging_team->members.count,
                                clear_team->members.list,
                                clear_team->members.count*sizeof(clear_team->members.list[0]));
                            nmerging_team->members.count += clear_team->members.count;
                            for (size_t x = 0; x < mt.members.count; ++x){
                                s_ENV.member_id_2_team_id[mt.members.list[x]] = mtid;
                            }
                            //free team, update team id ?
                            s_ENV.matching_team_pool.free(clear_team_id);
                        }
                        clear_bucket.waitings.list[merged_bucket_idxs[x].num_idx].\
                            teams.lremove(merged_bucket_idxs[x].pos_idx);
                    }
                    //insert as a total team (recursively merging)
                    merge_bucket.matchers.binsert(nmqst);
                }
            }
        }
    }
    return 0;
}

MatchingTeam_ST *    MatchingQueue::exit(uint64_t member_id){
    if (MATCHED_RESULT_TEAM_NUM > MQ.results.count){
        this->update(0, 500);//update for matched result fastly
    }
    size_t team_id = find_team_id_by_member_id(member_id);
    if (!team_id){ //not found
        return NULL;
    }
    MatchingTeam_ST * ptr_team = s_ENV.matching_team_pool.ptr(team_id);
    if (!ptr_team){
        return NULL;
    }
    MatchingBucket_ST & bucket = MQ.buckets[ptr_team->point.level];
    MatchingQueueMatcher_ST mqms;
    mqms.construct();
    mqms.team_elo = ptr_team->point.elo;
    mqms.team_id = team_id;
    int fid =  bucket.waitings.list[ptr_team->members.count-1].teams.lfind(mqms);
    if (fid >= 0){
        for (int i = fid; i < bucket.waitings.list[ptr_team->members.count - 1].teams.count; ++i){
            if (bucket.waitings.list[ptr_team->members.count - 1].teams.list[i].team_id == team_id){
                bucket.waitings.list[ptr_team->members.count - 1].teams.lremove(i);
                return ptr_team;
            }
        }
    }
    fid = bucket.matchers.lfind(mqms);
    if (fid >= 0){
        for (int i = fid; i < bucket.matchers.count; ++i){
            if (bucket.matchers.list[i].team_id == team_id){
                bucket.matchers.lremove(i);
                return ptr_team;
            }
        }
    }
    //result no should exit
    return NULL;
}
int     MatchingQueue::update(int past_ms, int checked_num){
    MQ.cur_ms_insec += past_ms;
    if (MQ.cur_ms_insec >= 1000){
        MQ.cur_time += MQ.cur_ms_insec / 1000;
        MQ.cur_ms_insec %= 1000;
    }
    //matching
    int iResultAvail = MACTCHING_QUEUE_MAX_MACTCHED_RESULT_NUM - MQ.results.count;
    int iCheckedNum = 0;
    for (int i = 0; i <= MACTCHING_QUEUE_MAX_LEVEL; ++i){
        MatchingBucket_ST &    bucket = MQ.buckets[i];
        MatchedResult_ST *     result = &MQ.results.list[MQ.results.count];
        result->teams.count = MATCHED_RESULT_TEAM_NUM;
        for (int j = bucket.matchers.count - 1; j >= 1 && iResultAvail > 0; j -= 2){
            result->teams[MATCHED_RESULT_TEAM_L] = bucket.matchers[j];
            result->teams[MATCHED_RESULT_TEAM_R] = bucket.matchers[j - 1];
            printf("same level:%d num:%d matched : j=%d j+1=%d [%d] <-> [%d]\n", i,
                bucket.matchers.count, j, j - 1,
                bucket.matchers[j].team_elo,
                bucket.matchers[j - 1].team_elo);

            ++MQ.results.count;
            --iResultAvail;
            ++iCheckedNum;
            bucket.matchers.count -= 2;
        }
    }
    if (iCheckedNum >= checked_num){
        return iCheckedNum;
    }
    //check rest span a section
    bool bMatching = false; //matching state
    int  iLeftLv = -1;
    MatchingQueueMatcher_ST result_teams[MATCHED_RESULT_TEAM_NUM];
    memset(&result_teams, 0, sizeof(result_teams));
    for (int i = MACTCHING_QUEUE_MAX_LEVEL; i >= 0; --i){
        MatchingBucket_ST &    bucket = MQ.buckets[i];
        MatchedResult_ST *     result = &MQ.results.list[MQ.results.count];
        if (bucket.matchers.count > 0 && iResultAvail > 0){
            if (bMatching){ //
                //checking
                result_teams[MATCHED_RESULT_TEAM_R] = bucket.matchers.list[0];
                if (_check_result_matched(MQ.cur_time, result_teams)){
                    MQ.results.lappend();
                    result->teams[MATCHED_RESULT_TEAM_L] = result_teams[MATCHED_RESULT_TEAM_L].team_id;
                    result->teams[MATCHED_RESULT_TEAM_R] = result_teams[MATCHED_RESULT_TEAM_R].team_id;
                    ++MQ.results.count;
                    --iResultAvail;
                    bMatching = false;
                    //remove
                    assert(iLeftLv >= 0);
                    assert(MQ.buckets.list[iLeftLv].matchers.count == 1);
                    --MQ.buckets.list[iLeftLv].matchers.count;

                    assert(bucket.matchers.count == 1);
                    --bucket.matchers.count;
                    printf("found span pair ok , then find next pair ...");
                    //find next left
                    i = iLeftLv;
                    iLeftLv = -1;
                }
            }
            else {
                result_teams[MATCHED_RESULT_TEAM_L] = bucket.matchers.list[0];
                iLeftLv = i;
                bMatching = true;
            }
        }
    }
    return iCheckedNum;
}
const MatchedResult_ST * MatchingQueue::fetch_results(int * result_num){

}

