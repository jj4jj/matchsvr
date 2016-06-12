#include <time.h>
#include <unordered_map>
#include <algorithm>
#include <cassert>
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

struct MatchingPoolImpl {
    MatchingQueue_ST queue;
};
MatchingPool::MatchingPool(){
    impl_ = new MatchingPoolImpl();
}

MatchingPool::~MatchingPool(){
    if (impl_){
        delete impl_;
        impl_ = NULL;
    }
}
int     MatchingPool::update_matching_object(uint64_t id, int elo, int lv){
    MatchingObject_ST mo;
    mo.construct();
    mo.id = id;
    mo.point.elo = elo;
    mo.point.level = lv;
    MatchingObject_ST * pf = s_ENV.matching_object_pool.find(mo);
    if (pf){
        pf->point = mo.point;
        return 0;
    }
    if (s_ENV.matching_object_pool.insert(mo)){
        return 0;
    }
    return -1;
}
static inline size_t find_team_id_by_member_id(uint64_t member_id){
    auto it = s_ENV.member_id_2_team_id.find(member_id);
    if (it != s_ENV.member_id_2_team_id.end()){
        return it->second;
    }
    return 0;
}
#define MQ  (impl_->queue)
int     MatchingPool::init(int team_member_num){
    matching_service_init();
    //////////////////////////////////////////////
    impl_->queue.construct();
    impl_->queue.team_member_num = team_member_num;
    impl_->queue.cur_time = ::time(NULL);
    MQ.buckets.count = MACTCHING_QUEUE_MAX_LEVEL_NUM;
    for (int i = 0; i < MACTCHING_QUEUE_MAX_LEVEL_NUM; ++i){
        MQ.buckets[i].waitings.count = MATCHING_TEAM_MAX_WAITING_MERGING_MEMBER_NUM;
        for (size_t j = 0; j < MQ.buckets[i].waitings.count; ++j){
            MQ.buckets[i].waitings[j].teams.construct();
        }
    }
    return 0;
}
static inline void _elo_team_init(MatchingTeam_ST & t){
    int point_add_extra_percent = 0;
    switch (t.members.count){
    case 1:
        break;
    case 2:
        point_add_extra_percent = 5;
        break;
    case 3:
        point_add_extra_percent = 10;
        break;
    default:
        break;
    }
    MatchingObject_ST k;
    size_t elo_product = 1;
    size_t lv_product = 1;
    for (size_t i = 0; i < t.members.count; ++i){
        k.id = t.members[i];
        MatchingObject_ST * pObj = s_ENV.matching_object_pool.find(k);
        assert(pObj);
        elo_product *= pObj->point.elo;
        lv_product *= pObj->point.level;
    }
    if (t.members.count > 1){
        t.point.elo = pow((float)elo_product, 1.0 / (float)(t.members.count));
        t.point.level = pow((float)lv_product, 1.0 / (float)(t.members.count));
    }
    else {
        t.point.elo = elo_product;
        t.point.level = lv_product;
    }
    if (point_add_extra_percent > 0){
        t.point.elo *= (100 + point_add_extra_percent);
        t.point.elo /= 100;
    }
}
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
static inline bool _check_team_waiting_can_merged(uint32_t cur_time, const MatchingTeam_ST & t1, const MatchingQueueMatcher_ST & t2){
    MatchingTeam_ST * rt = s_ENV.matching_team_pool.ptr(t2.team_id);
    assert(rt);
    int wait_time = std::min(cur_time - t1.join_time, cur_time - rt->join_time);
    int lv_diff = t1.point.level - rt->point.level;
    if (lv_diff < 0){ lv_diff = -lv_diff; }
    if (wait_time >= lv_diff * 15){
        return true;
    }
    return false;
}
static inline bool _check_result_matched(uint32_t cur_time, MatchingQueueMatcher_ST results[MATCHED_RESULT_TEAM_NUM]){
    MatchingTeam_ST * lt = s_ENV.matching_team_pool.ptr(results[MATCHED_RESULT_TEAM_L].team_id);
    MatchingTeam_ST * rt = s_ENV.matching_team_pool.ptr(results[MATCHED_RESULT_TEAM_R].team_id);
    assert(lt && rt);
    int llv = lt->point.level;// _elo_to_level(results[MATCHED_RESULT_TEAM_L].team_elo);
    int rlv = rt->point.level;// _elo_to_level(results[MATCHED_RESULT_TEAM_R].team_elo);
    int wait_time = std::min(cur_time - lt->join_time, cur_time - rt->join_time);
    int lv_diff = llv - rlv;
    if (lv_diff < 0){ lv_diff = -lv_diff; }
    if (wait_time >= lv_diff * 15){
        printf("check result matched success -> lv diff:%d wait min time:%d [%u: %u %u in time %u]<->[%u: %u %u in time %u]\n",
            lv_diff, wait_time,
            results[MATCHED_RESULT_TEAM_L].team_id, lt->point.elo, llv, cur_time - lt->join_time,
            results[MATCHED_RESULT_TEAM_R].team_id, rt->point.elo, rlv, cur_time - rt->join_time);
        return true;
    }
    return false;
}
static inline void _check_merge_team(MatchingPoolImpl * impl_){
    MatchingTeam_ST merged;
    merged.construct();
    size_t merged_buckets_num = 0;
    struct {
        int level;
        int num_idx;
        int pos_idx;
    } merged_bucket_idxs[MATCHING_TEAM_MAX_MEMBER_NUM];
    //check merging all
    for (int level = MACTCHING_QUEUE_MAX_LEVEL; level >= 0; --level){
        MatchingBucket_ST &    mqbucket = MQ.buckets[level];
        for (int wq_nmem = MQ.team_member_num - 1; wq_nmem > 0; --wq_nmem){ //member number max -> min [N-1 => 1]
            int nmem_idx = wq_nmem - 1;
            MatchingWaitigTeam_ST & wmt = mqbucket.waitings[nmem_idx];
            for (int k = wmt.teams.count - 1; k >= 0; --k){
                //team can join  check ?
                if (merged_buckets_num > 0 && !_check_team_waiting_can_merged(MQ.cur_time, merged, wmt.teams[k])){
                    continue;
                }
                if (merged.members.count + wq_nmem <= MQ.team_member_num){
                    //merge team
                    merged_bucket_idxs[merged_buckets_num].level = level;
                    merged_bucket_idxs[merged_buckets_num].num_idx = nmem_idx;
                    merged_bucket_idxs[merged_buckets_num].pos_idx = k;
                    ++merged_buckets_num;
                    ///////////////////////////////////////////////////////
                    if (merged.members.count == 0){
                        merged.point.elo = wmt.teams[k].team_elo;
                        merged.point.level = level;
                        merged.members.count = wq_nmem;
                    }
                    else {
                        merged.point.elo = sqrt(merged.point.elo * wmt.teams[k].team_elo);
                        merged.point.level = sqrt(merged.point.level * level);
                        merged.members.count += wq_nmem;
                    }
                }
                //found one
                if (merged.members.count == MQ.team_member_num){
                    printf("found one merge matching team with %zu teams [%u]\n",
                        merged_buckets_num, merged.point.elo);
                    MatchingBucket_ST & merge_bucket = MQ.buckets[merged.point.level];
                    if (merge_bucket.matchers.count >= MATCHING_TEAM_MAX_NUM_IN_LEVEL){
                        memset(&merged, 0, sizeof(merged));
                        merged_buckets_num = 0;
                        continue; //matching queueu is full
                    }
                    //new matcher
                    MatchingQueueMatcher_ST nmqst;
                    nmqst.construct();
                    nmqst.team_elo = merged.point.elo;
                    MatchingTeam_ST * nmerging_team = NULL;
                    printf("clear merged %zu teams\n", merged_buckets_num);
                    //clear merged
                    for (size_t x = 0; x < merged_buckets_num; ++x){
                        MatchingBucket_ST & clear_bucket = MQ.buckets[merged_bucket_idxs[x].level];
                        MatchingWaitigTeam_ST & clear_wq = clear_bucket.waitings[merged_bucket_idxs[x].num_idx];
                        size_t clear_team_id = clear_wq.teams[merged_bucket_idxs[x].pos_idx].team_id;
                        MatchingTeam_ST *  clear_team = s_ENV.matching_team_pool.ptr(clear_team_id);
                        assert(clear_team);
                        if (x == 0){
                            nmqst.team_id = clear_team_id;
                            nmerging_team = clear_team;
                            nmerging_team->join_time = MQ.cur_time;
                            nmerging_team->point = merged.point;
                        }
                        else {
                            memcpy(nmerging_team->members.list + nmerging_team->members.count,
                                clear_team->members.list,
                                clear_team->members.count*sizeof(clear_team->members[0]));
                            nmerging_team->members.count += clear_team->members.count;
                            for (size_t x = 0; x < clear_team->members.count; ++x){ //update team map
                                s_ENV.member_id_2_team_id[clear_team->members[x]] = nmqst.team_id;
                            }
                            //free team, update team id ?
                            printf("free no:%d team:%zu\n", __LINE__, clear_team_id);
                            s_ENV.matching_team_pool.free(clear_team_id);
                        }
                        //here no idx change when removing for pos idx is reversed
                        //when in same nidx (member num) , it's safe
                        printf("before remove waiting queue x=%zu [wq count=%u]\n", x,
                            clear_wq.teams.count);
                        clear_wq.teams.lremove(merged_bucket_idxs[x].pos_idx);
                        printf("after remove waiting queue x=%zu "
                            "level=%u num=%u pos=%u [wq count=%u]\n", x,
                            merged_bucket_idxs[x].level,
                            merged_bucket_idxs[x].num_idx + 1,
                            merged_bucket_idxs[x].pos_idx,
                            clear_wq.teams.count);
                    }
                    //clear merged temp
                    memset(&merged, 0, sizeof(merged));
                    merged_buckets_num = 0;
                    //insert as a total team (recursively merging)
                    merge_bucket.matchers.binsert(nmqst);
                }
            }
        }
    }
}
int     MatchingPool::join(const std::vector<uint64_t> & members){
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
    assert(s_ENV.matching_team_pool.ptr(mtid));
    *s_ENV.matching_team_pool.ptr(mtid) = mt;

    int matching_lv = mt.point.level;
    MatchingBucket_ST & bucket = MQ.buckets[matching_lv];

    MatchingQueueMatcher_ST mqst;
    mqst.construct();
    mqst.team_elo = mt.point.elo;
    mqst.team_id = mtid;

    printf("join matching team %zu#%u [%u %u] (", mtid,  mt.members.count, matching_lv, mt.point.elo);
    for (size_t x = 0; x < mt.members.count; ++x){
        printf("%lu,", mt.members[x]);
    }
    printf(")\n");


    if (mt.members.count == MQ.team_member_num){
        //insert into matchers
        if (bucket.matchers.count >= MATCHING_TEAM_MAX_NUM_IN_LEVEL){
            printf("free no:%d team:%zu\n", __LINE__, mtid);
            s_ENV.matching_team_pool.free(mtid);
            return -3;
        }
        if (bucket.matchers.binsert(mqst)){
            printf("free no:%d team:%zu\n", __LINE__, mtid);
            s_ENV.matching_team_pool.free(mtid);
            return -4;
        }
        for (size_t x = 0; x < mt.members.count; ++x){
            s_ENV.member_id_2_team_id[mt.members[x]] = mtid;
        }
    }
    else {//waiting merging
        if (bucket.waitings[mt.members.count - 1].teams.count == MACTCHING_BUCKET_MAX_WAITING_TEAM_NUM){
            printf("free no:%d team:%zu\n", __LINE__, mtid);
            s_ENV.matching_team_pool.free(mtid);
            return -5; //waiting queue is also full
        }
        if (bucket.waitings[mt.members.count - 1].teams.binsert(mqst)){
            printf("free no:%d team:%zu\n", __LINE__, mtid);
            s_ENV.matching_team_pool.free(mtid);
            return -6;
        }
        for (size_t x = 0; x < mt.members.count; ++x){
            s_ENV.member_id_2_team_id[mt.members[x]] = mtid;
        }
    }
    _check_merge_team(impl_);
    return 0;
}
static inline void _free_team(size_t team_id){
    MatchingTeam_ST * team = s_ENV.matching_team_pool.ptr(team_id);
    printf("free no:%d team:%zu ptr:%p\n", __LINE__, team_id, team);
    for (size_t k = 0; k < team->members.count; ++k){
        s_ENV.member_id_2_team_id.erase(team->members[k]);
    }
    s_ENV.matching_team_pool.free(team_id);
}
const MatchingTeam_ST *  MatchingPool::get_team(size_t team_id){
    return s_ENV.matching_team_pool.ptr(team_id);
}
void                 MatchingPool::free_team(const MatchingTeam_ST * team){
    size_t team_id = s_ENV.matching_team_pool.id((MatchingTeam_ST *)team);
    printf("free no:%d team:%zu\n", __LINE__, team_id);
    assert(team_id > 0);
    _free_team(team_id);
}
const MatchingTeam_ST *    MatchingPool::quit(uint64_t member_id){
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
    if (ptr_team->members.count < MQ.team_member_num){
        int inm = ptr_team->members.count - 1;
        int fid = bucket.waitings[inm].teams.bfind(mqms);
        if (fid >= 0){
            for (size_t i = fid; i < bucket.waitings[inm].teams.count; ++i){
                if (bucket.waitings[inm].teams[i].team_id == team_id){
                    bucket.waitings[inm].teams.lremove(i);
                    return ptr_team;
                }
            }
        }
        assert(false);
    }
    else {
        int fid = bucket.matchers.bfind(mqms);
        if (fid >= 0){
            for (size_t i = fid; i < bucket.matchers.count; ++i){
                if (bucket.matchers[i].team_id == team_id){
                    bucket.matchers.lremove(i);
                    return ptr_team;
                }
            }
        }
        //should be in result queue
        return NULL;
    }
}
uint32_t MatchingPool::time() const {
    return MQ.cur_time;
}
int     MatchingPool::update(int past_ms, int max_checked_num){
    MQ.cur_ms_insec += past_ms;
    if (MQ.cur_ms_insec >= 1000){
        MQ.cur_time += MQ.cur_ms_insec / 1000;
        MQ.cur_ms_insec %= 1000;
    }
    //matching
    int result_avail = MACTCHING_QUEUE_MAX_MACTCHED_RESULT_NUM - MQ.results.count;
    int checked_num = 0;
    MatchedResult_ST result;
    result.construct();
    result.teams.count = MATCHED_RESULT_TEAM_NUM;
    for (int i = 0; i <= MACTCHING_QUEUE_MAX_LEVEL; ++i){
        MatchingBucket_ST &    bucket = MQ.buckets[i];
        for (int j = bucket.matchers.count - 1; j >= 1 && result_avail > 0; j -= 2){
            result.teams[MATCHED_RESULT_TEAM_L] = bucket.matchers[j].team_id;
            result.teams[MATCHED_RESULT_TEAM_R] = bucket.matchers[j - 1].team_id;

            printf("same level:%d num:%u matched : [%u:%u] <-> [%u:%u]\n", i,
                bucket.matchers.count,
                bucket.matchers[j].team_id, bucket.matchers[j].team_elo,
                bucket.matchers[j - 1].team_id, bucket.matchers[j - 1].team_elo);

            MQ.results.lappend(result);
            --result_avail;
            ++checked_num;
            bucket.matchers.count -= 2;
        }
    }
    if (checked_num >= max_checked_num){
        return checked_num;
    }
    //check rest span a section
    bool matching = false; //matching state
    int  left_level = -1;
    MatchingQueueMatcher_ST result_teams[MATCHED_RESULT_TEAM_NUM];
    memset(&result_teams, 0, sizeof(result_teams));
    result.teams[MATCHED_RESULT_TEAM_L] = 0;
    result.teams[MATCHED_RESULT_TEAM_R] = 0;
    for (int i = MACTCHING_QUEUE_MAX_LEVEL; i >= 0; --i){
        MatchingBucket_ST &    bucket = MQ.buckets[i];
        if (bucket.matchers.count > 0 && result_avail > 0 && checked_num < max_checked_num){
            assert(bucket.matchers[0].team_id);
            if (matching){ //
                //checking
                result_teams[MATCHED_RESULT_TEAM_R] = bucket.matchers[0];
                if (_check_result_matched(MQ.cur_time, result_teams)){
                    result.teams[MATCHED_RESULT_TEAM_L] = result_teams[MATCHED_RESULT_TEAM_L].team_id;
                    result.teams[MATCHED_RESULT_TEAM_R] = result_teams[MATCHED_RESULT_TEAM_R].team_id;
                    MQ.results.lappend(result);
                    --result_avail;
                    matching = false;
                    //
                    printf("span level diff:%d matched : [%u:%u] <-> [%u:%u]\n", i - left_level,
                        result_teams[MATCHED_RESULT_TEAM_L].team_id, result_teams[MATCHED_RESULT_TEAM_L].team_elo,
                        result_teams[MATCHED_RESULT_TEAM_R].team_id, result_teams[MATCHED_RESULT_TEAM_R].team_elo);
                    
                    //clear team result
                    result.teams[MATCHED_RESULT_TEAM_L] = 0;
                    result.teams[MATCHED_RESULT_TEAM_R] = 0;
                    memset(&result_teams, 0, sizeof(result_teams));

                    //remove index
                    assert(left_level >= 0);
                    assert(MQ.buckets[left_level].matchers.count == 1);
                    --MQ.buckets[left_level].matchers.count;

                    assert(bucket.matchers.count == 1);
                    --bucket.matchers.count;
                    printf("found span pair ok , then find next pair ...\n");
                    //find next left
                    i = left_level;
                    left_level = -1;
                }
            }
            else {
                result_teams[MATCHED_RESULT_TEAM_L] = bucket.matchers[0];
                left_level = i;
                matching = true;
            }
        }
    }
    return checked_num;
}
const MatchedResult_ST * MatchingPool::fetch_results(size_t * result_num){
    assert(result_num );
    if (*result_num == 0){
        return NULL;
    }
    //*result_num = MQ.results.count;
    assert(MQ.results.count >= MQ.fetched_results);
    MatchedResult_ST * p = MQ.results.list + MQ.fetched_results;
    if (*result_num > (MQ.results.count - MQ.fetched_results)){//cut
        *result_num = MQ.results.count - MQ.fetched_results;
    }
    MQ.fetched_results += *result_num;
    return p;
}
void                     MatchingPool::clear_fetched_results(){
    if (MQ.fetched_results == 0){
        return;
    }
    assert(MQ.results.count >= MQ.fetched_results);
    for (uint32_t i = 0; i < MQ.fetched_results; ++i){
        for (int j = 0; j < MATCHED_RESULT_TEAM_NUM; ++j){
            _free_team(MQ.results[i].teams[j]);
        }
    }
    //
    if (MQ.results.count > MQ.fetched_results){
        memmove(MQ.results.list, MQ.results.list + MQ.fetched_results,
            (MQ.results.count - MQ.fetched_results)*sizeof(MQ.results.list[0]));
    }
    MQ.results.count -= MQ.fetched_results;
    printf("clear result num:%u rest:%u\n", MQ.fetched_results, MQ.results.count);
    MQ.fetched_results = 0;
}

