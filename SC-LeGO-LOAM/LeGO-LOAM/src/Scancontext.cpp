#include "Scancontext.h"

// namespace SC2
// {

void coreImportTest (void)
{
    cout << "scancontext lib is successfully imported." << endl;
} // coreImportTest


float rad2deg(float radians)
{
    return radians * 180.0 / M_PI;
}

float deg2rad(float degrees)
{
    return degrees * M_PI / 180.0;
}


float xy2theta( const float & _x, const float & _y )
{
    if ( _x >= 0 & _y >= 0) 
        return (180/M_PI) * atan(_y / _x);

    if ( _x < 0 & _y >= 0) 
        return 180 - ( (180/M_PI) * atan(_y / (-_x)) );

    if ( _x < 0 & _y < 0) 
        return 180 + ( (180/M_PI) * atan(_y / _x) );

    if ( _x >= 0 & _y < 0)
        return 360 - ( (180/M_PI) * atan((-_y) / _x) );
} // xy2theta


MatrixXd circshift( MatrixXd &_mat, int _num_shift )
{
    // shift columns to right direction 
    assert(_num_shift >= 0);

    if( _num_shift == 0 )
    {
        MatrixXd shifted_mat( _mat );
        return shifted_mat; // Early return 
    }

    MatrixXd shifted_mat = MatrixXd::Zero( _mat.rows(), _mat.cols() );
    for ( int col_idx = 0; col_idx < _mat.cols(); col_idx++ )
    {
        int new_location = (col_idx + _num_shift) % _mat.cols();
        shifted_mat.col(new_location) = _mat.col(col_idx);
    }

    return shifted_mat;

} // circshift


std::vector<float> eig2stdvec( MatrixXd _eigmat )
{
    std::vector<float> vec( _eigmat.data(), _eigmat.data() + _eigmat.size() );
    return vec;
} // eig2stdvec

// 回环检测第二步精匹配 中的 第2小步精确匹配（计算当前偏移量下的距离 = 1 - 两个向量的cos值 （两向量越相似，越重合，cos值越大，此距离越小））
double SCManager::distDirectSC ( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    int num_eff_cols = 0; // i.e., to exclude all-nonzero sector
    double sum_sector_similarity = 0;
    for ( int col_idx = 0; col_idx < _sc1.cols(); col_idx++ )
    {
        VectorXd col_sc1 = _sc1.col(col_idx);
        VectorXd col_sc2 = _sc2.col(col_idx);
        
        if( col_sc1.norm() == 0 | col_sc2.norm() == 0 )
            continue; // don't count this sector pair. 

        double sector_similarity = col_sc1.dot(col_sc2) / (col_sc1.norm() * col_sc2.norm());

        sum_sector_similarity = sum_sector_similarity + sector_similarity;
        num_eff_cols = num_eff_cols + 1;
    }
    
    double sc_sim = sum_sector_similarity / num_eff_cols;
    return 1.0 - sc_sim;

} // distDirectSC

// 回环检测第二步精匹配 中的 第1小步快速的粗对齐（用Vector Key（存储着每列均值的向量）的norm作差，差越小说明向量越相似，实现快速对齐）
int SCManager::fastAlignUsingVkey( MatrixXd & _vkey1, MatrixXd & _vkey2)
{
    int argmin_vkey_shift = 0;
    double min_veky_diff_norm = 10000000;
    for ( int shift_idx = 0; shift_idx < _vkey1.cols(); shift_idx++ )
    {
        MatrixXd vkey2_shifted = circshift(_vkey2, shift_idx);

        MatrixXd vkey_diff = _vkey1 - vkey2_shifted;

        double cur_diff_norm = vkey_diff.norm();
        if( cur_diff_norm < min_veky_diff_norm )
        {
            argmin_vkey_shift = shift_idx;
            min_veky_diff_norm = cur_diff_norm;
        }
    }

    return argmin_vkey_shift;

} // fastAlignUsingVkey

// 找到query和target之间的最优偏移量
std::pair<double, int> SCManager::distanceBtnScanContext( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    // 1. fast align using variant key (not in original IROS18)  第一步，快速的粗对齐
    MatrixXd vkey_sc1 = makeSectorkeyFromScancontext( _sc1 );
    MatrixXd vkey_sc2 = makeSectorkeyFromScancontext( _sc2 );
    int argmin_vkey_shift = fastAlignUsingVkey( vkey_sc1, vkey_sc2 );   // 用Vector Key（存储着每列均值的向量）的norm作差，差越小说明向量越相似，实现快速对齐

    const int SEARCH_RADIUS = round( 0.5 * SEARCH_RATIO * _sc1.cols() ); // a half of search range 
    std::vector<int> shift_idx_search_space { argmin_vkey_shift };
    for ( int ii = 1; ii < SEARCH_RADIUS + 1; ii++ )    // 毕竟上面是快速搜索，可能不准，将取得的偏移量的左右几个作为小搜索范围，用于接下来的精细搜索
    {
        shift_idx_search_space.push_back( (argmin_vkey_shift + ii + _sc1.cols()) % _sc1.cols() );
        shift_idx_search_space.push_back( (argmin_vkey_shift - ii + _sc1.cols()) % _sc1.cols() );
    }
    std::sort(shift_idx_search_space.begin(), shift_idx_search_space.end());

    // 2. fast columnwise diff 第二步，精细对齐
    int argmin_shift = 0;
    double min_sc_dist = 10000000;
    for ( int num_shift: shift_idx_search_space )
    {
        MatrixXd sc2_shifted = circshift(_sc2, num_shift);          // 获得_sc2偏移num_shift后的结果
        double cur_sc_dist = distDirectSC( _sc1, sc2_shifted );     // 计算当前偏移量下的距离 = 1 - 两个向量的cos值 （两向量越相似，越重合，cos值越大，此距离越小）
        if( cur_sc_dist < min_sc_dist )                             // 找到最小的距离，更新此时的偏移量
        {
            argmin_shift = num_shift;
            min_sc_dist = cur_sc_dist;
        }
    }

    return make_pair(min_sc_dist, argmin_shift);     // 返回最小距离和对应的偏移量

} // distanceBtnScanContext

// 创建scan context
MatrixXd SCManager::makeScancontext( pcl::PointCloud<SCPointType> & _scan_down )
{
    TicToc t_making_desc;

    int num_pts_scan_down = _scan_down.points.size();

    // main
    const int NO_POINT = -1000;
    MatrixXd desc = NO_POINT * MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);

    SCPointType pt;
    float azim_angle, azim_range; // wihtin 2d plane
    int ring_idx, sctor_idx;
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z = _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // naive adding is ok (all points should be > 0).

        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);

        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )
            continue;
        // 计算ring index和 sector index
        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );

        // taking maximum z 
        if ( desc(ring_idx-1, sctor_idx-1) < pt.z ) // -1 means cpp starts from 0
            desc(ring_idx-1, sctor_idx-1) = pt.z; // update for taking maximum value at that bin
    }

    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )
            if( desc(row_idx, col_idx) == NO_POINT )
                desc(row_idx, col_idx) = 0;

    t_making_desc.toc("PolarContext making");

    return desc;
} // SCManager::makeScancontext

// 计算ring key(取每一行的均值)，用于后续回环检测第一步构建kd-tree并进行粗匹配
MatrixXd SCManager::makeRingkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: rowwise mean vector
    */
    Eigen::MatrixXd invariant_key(_desc.rows(), 1);
    for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ )
    {
        Eigen::MatrixXd curr_row = _desc.row(row_idx);
        invariant_key(row_idx, 0) = curr_row.mean();
    }

    return invariant_key;
} // SCManager::makeRingkeyFromScancontext

// 计算sector key(取每一行的均值),用于后续回环检测第二步精匹配中的粗对齐
MatrixXd SCManager::makeSectorkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: columnwise mean vector
    */
    Eigen::MatrixXd variant_key(1, _desc.cols());
    for ( int col_idx = 0; col_idx < _desc.cols(); col_idx++ )
    {
        Eigen::MatrixXd curr_col = _desc.col(col_idx);
        variant_key(0, col_idx) = curr_col.mean();
    }

    return variant_key;
} // SCManager::makeSectorkeyFromScancontext


void SCManager::makeAndSaveScancontextAndKeys( pcl::PointCloud<SCPointType> & _scan_down )
{
    Eigen::MatrixXd sc = makeScancontext(_scan_down); // v1 
    Eigen::MatrixXd ringkey = makeRingkeyFromScancontext( sc );
    Eigen::MatrixXd sectorkey = makeSectorkeyFromScancontext( sc );
    std::vector<float> polarcontext_invkey_vec = eig2stdvec( ringkey );

    polarcontexts_.push_back( sc ); 
    polarcontext_invkeys_.push_back( ringkey );
    polarcontext_vkeys_.push_back( sectorkey );
    polarcontext_invkeys_mat_.push_back( polarcontext_invkey_vec );

    // cout <<polarcontext_vkeys_.size() << endl;

} // SCManager::makeAndSaveScancontextAndKeys

// 回环检测的主函数
// 1. 计算当前scan context的ringkey（行均值向量）
// 2. 在ringkey树中找到最近的几个候选scan context（根据ringkey的距离）
// 3. 对每个候选scan context，计算其与当前scan context的距离（sectorkey的粗匹配 + cos距离的精匹配）（该距离值越小，表示越相似）
// 4. 如果找到的候选scan context的距离小于某个阈值，则认为存在loop，返回该候选scan context的ID和yaw的偏移角度

// ┌─────────────────────────────────────┐
// │ 获取当前帧 ringkey 和 描述符           │
// └──────────┬──────────────────────────┘
//            ↓
// ┌─────────────────────────────────────┐
// │ Step 1: KD-树搜索（粗筛选）           │
// │ 输入：当前 ringkey                   │
// │ 输出：10 个最近邻候选帧ID              │
// └──────────┬──────────────────────────┘
//            ↓
// ┌─────────────────────────────────────┐
// │ Step 2: 逐个精细对齐（精匹配）          │
// │ 对每个候选帧：                        │
// │  - 计算完整扫描上下文距离              │
// │  - 找最优yaw对齐角度                  │
// │  - 记录最小距离及对应帧ID              │
// └──────────┬──────────────────────────┘
//            ↓
// ┌─────────────────────────────────────┐
// │ Step 3: 阈值判断                      │
// │ 如果 min_dist < SC_DIST_THRES        │
// │   → 返回 (帧ID, yaw偏移)              │
// │ 否则 → 返回 (-1, 0.0)                │
// └─────────────────────────────────────┘
std::pair<int, float> SCManager::detectLoopClosureID ( void )
{
    // 初始化与获取当前帧数据
    int loop_id { -1 }; // init with -1, -1 means no loop (== LeGO-LOAM's variable "closestHistoryFrameID")

    auto curr_key = polarcontext_invkeys_mat_.back(); // current observation (query)
    auto curr_desc = polarcontexts_.back(); // current observation (query)

    /* 
     * step 1: candidates from ringkey tree_  
        第一步 粗筛选：用ringkey从KD tree中找到几个最近邻候选帧
     */
    // 检查数据量是否足够:排除最近的N帧（因为最近的NUM_EXCLUDE_RECENT帧不太可能是loop，且会占用树的搜索资源），如果不够则直接返回没有loop
    if( polarcontext_invkeys_mat_.size() < NUM_EXCLUDE_RECENT + 1)
    {
        std::pair<int, float> result {loop_id, 0.0};
        return result; // Early return 
    }

    // tree_ reconstruction (not mandatory to make everytime) 周期性用ringkey构建KD数：每隔TREE_MAKING_PERIOD_帧重新构建KD树，节省计算成本
    if( tree_making_period_conter % TREE_MAKING_PERIOD_ == 0) // to save computation cost
    {
        TicToc t_tree_construction;
        // 假设历史数据：
        //polarcontext_invkeys_mat_ = [key_0, key_1, ..., key_98, key_99]  // 共100帧
        // NUM_EXCLUDE_RECENT = 20

        // 分配后（排除最近的20帧）：
        // polarcontext_invkeys_to_search_ = [key_0, key_1, ..., key_79]  // 只有80帧
        polarcontext_invkeys_to_search_.clear();
        polarcontext_invkeys_to_search_.assign( polarcontext_invkeys_mat_.begin(), polarcontext_invkeys_mat_.end() - NUM_EXCLUDE_RECENT ) ;

        polarcontext_tree_.reset(); 
        polarcontext_tree_ = std::make_unique<InvKeyTree>(PC_NUM_RING /* dim */, polarcontext_invkeys_to_search_, 10 /* max leaf */ );
        // tree_ptr_->index->buildIndex(); // inernally called in the constructor of InvKeyTree (for detail, refer the nanoflann and KDtreeVectorOfVectorsAdaptor)
        t_tree_construction.toc("Tree construction");
    }
    tree_making_period_conter = tree_making_period_conter + 1;
        
    double min_dist = 10000000; // init with somthing large
    int nn_align = 0;
    int nn_idx = 0;

    // knn search
    std::vector<size_t> candidate_indexes( NUM_CANDIDATES_FROM_TREE ); 
    std::vector<float> out_dists_sqr( NUM_CANDIDATES_FROM_TREE );

    TicToc t_tree_search;
    nanoflann::KNNResultSet<float> knnsearch_result( NUM_CANDIDATES_FROM_TREE );    // 创建一个容器，存放10个近邻结果
    knnsearch_result.init( &candidate_indexes[0], &out_dists_sqr[0] );              // 初始化结果容器的指针
    polarcontext_tree_->index->findNeighbors( knnsearch_result, &curr_key[0] /* query */, nanoflann::SearchParams(10) ); // 在树中搜索，自动填充 candidate_indexes 和 out_dists_sqr
    t_tree_search.toc("Tree search");

    /* 
     *  step 2: pairwise distance (find optimal columnwise best-fit using cosine distance) 
        第二步 精匹配：经过第一步找到最相似的候选帧后，将query与候选帧逐个进行精匹配，找到最小距离对应的候选帧ID和yaw偏移量（即最佳对齐角度）。注意：这里的距离是scan context之间的距离，而不是ringkey之间的距离。
     */
    TicToc t_calc_dist;   
    for ( int candidate_iter_idx = 0; candidate_iter_idx < NUM_CANDIDATES_FROM_TREE; candidate_iter_idx++ )
    {
        MatrixXd polarcontext_candidate = polarcontexts_[ candidate_indexes[candidate_iter_idx] ];  // 取出候选帧
        std::pair<double, int> sc_dist_result = distanceBtnScanContext( curr_desc, polarcontext_candidate ); 
        
        double candidate_dist = sc_dist_result.first;
        int candidate_align = sc_dist_result.second;

        if( candidate_dist < min_dist )
        {
            min_dist = candidate_dist;
            nn_align = candidate_align;

            nn_idx = candidate_indexes[candidate_iter_idx];
        }
    }
    t_calc_dist.toc("Distance calc");

    /* 
     * loop threshold check
     */
    if( min_dist < SC_DIST_THRES )
    {
        loop_id = nn_idx; 
    
        // std::cout.precision(3);  例：第500帧与第120帧形成回环，需要旋转45°对齐。
        cout << "[Loop found] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        cout << "[Loop found] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }
    else
    {
        std::cout.precision(3); 
        cout << "[Not loop] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        cout << "[Not loop] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }

    // To do: return also nn_align (i.e., yaw diff)
    float yaw_diff_rad = deg2rad(nn_align * PC_UNIT_SECTORANGLE);
    std::pair<int, float> result {loop_id, yaw_diff_rad};

    return result;

} // SCManager::detectLoopClosureID

// } // namespace SC2