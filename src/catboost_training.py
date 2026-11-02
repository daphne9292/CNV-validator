import joblib
import numpy as np
import catboost as cb
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report
import pandas as pd
import os

f_list_drop1 = ['svtype', 'report_len','ave_normal_baf_dev_of_interval',\
                'ave_tumor_baf_dev_leftside_of_interval', 'ave_normal_baf_dev_leftside_of_interval', \
                'ave_tumor_baf_dev_rightside_of_interval', 'ave_normal_baf_dev_rightside_of_interval',\
               'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval',\
                "SR_SA_cluster_num_n_lbp","SR_SA_cluster_num_n_lbp", \
                'SA_in_rbp_neigb_divide_SR_n_lbp', 'SA_in_lbp_neigb_divide_SR_n_rbp', \
                'maxC_SRnum_divide_totalR_n_rbp', 'maxC_SRnum_divide_totalR_n_lbp', \
                'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',\
                "SR_SA_cluster_num_t_lbp","SR_SA_cluster_num_t_lbp", \
                'SA_in_rbp_neigb_divide_SR_t_lbp', 'SA_in_lbp_neigb_divide_SR_t_rbp', \
                'maxC_SRnum_divide_totalR_t_rbp', 'maxC_SRnum_divide_totalR_t_lbp', \
                'mate_in_rbp_neigb_divide_abRP_n_lbp', 'mate_in_lbp_neigb_divide_abRP_n_rbp', \
                'valid_t_RP_percent_left_bp_interval','valid_t_RP_percent_right_bp_interval',\
                'valid_n_RP_percent_left_bp_interval','valid_n_RP_percent_right_bp_interval', \
                'n_depth_ratio_of_leftside_with_interval', 'n_depth_ratio_of_rightside_with_interval', \
                'RP_mate_loci_cluster_num_t_lbp','RP_mate_loci_cluster_num_t_rbp',\
                'RP_mate_loci_cluster_num_n_lbp','RP_mate_loci_cluster_num_n_rbp']

f_list_drop2=['svtype', 'report_len','left_label', 'right_label','ave_normal_baf_dev_of_interval',\
              'ave_tumor_baf_dev_leftside_of_interval', 'ave_normal_baf_dev_leftside_of_interval', \
              'ave_tumor_baf_dev_rightside_of_interval', 'ave_normal_baf_dev_rightside_of_interval',\
              'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval',\
              "SR_SA_cluster_num_n_lbp","SR_SA_cluster_num_n_lbp", \
              'SA_in_rbp_neigb_divide_SR_n_lbp', 'SA_in_lbp_neigb_divide_SR_n_rbp', \
              'maxC_SRnum_divide_totalR_n_rbp','maxC_SRnum_divide_totalR_n_lbp', \
              'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval', \
              "SR_SA_cluster_num_t_lbp", "SR_SA_cluster_num_t_lbp", \
              'SA_in_rbp_neigb_divide_SR_t_lbp', 'SA_in_lbp_neigb_divide_SR_t_rbp', \
              'maxC_SRnum_divide_totalR_t_rbp', 'maxC_SRnum_divide_totalR_t_lbp', \
              'mate_in_rbp_neigb_divide_abRP_n_lbp','mate_in_lbp_neigb_divide_abRP_n_rbp', \
              'valid_t_RP_percent_left_bp_interval', 'valid_t_RP_percent_right_bp_interval', \
              'valid_n_RP_percent_left_bp_interval', 'valid_n_RP_percent_right_bp_interval', \
              'n_depth_ratio_of_leftside_with_interval', 'n_depth_ratio_of_rightside_with_interval', \
              'RP_mate_loci_cluster_num_t_lbp', 'RP_mate_loci_cluster_num_t_rbp', \
              'RP_mate_loci_cluster_num_n_lbp', 'RP_mate_loci_cluster_num_n_rbp']

class CNVCatBoost:
    def __init__(self):
        self.model = None

    def train(self, X, y):
        X_train, X_val, y_train, y_val = train_test_split(
            X, y,
            test_size=0.2,
            random_state=42,
            stratify=y
        )
        y_train = y_train.flatten()
        y_val = y_val.flatten()

        # ========== 样本类别权重（解决不均衡，和你原来逻辑一致） ==========
        cls_cnt = np.bincount(y_train.astype(int))
        total = len(y_train)
        cls_w = total / cls_cnt
        class_weights = cls_w.tolist()

        # 构建Pool
        train_pool = cb.Pool(data=X_train, label=y_train)
        val_pool = cb.Pool(data=X_val, label=y_val)
        print(class_weights)

        params = {
            "loss_function": "MultiClass",
            "eval_metric": "MultiClass",
            "boosting_type": "Plain",
            "learning_rate": 0.05,
            "depth": 6,
            "l2_leaf_reg": 0.9,
            "min_child_samples": 7,
            "class_weights": class_weights,   # 类别均衡权重
            "iterations": 200,
            "early_stopping_rounds": 50,
            "verbose": 50,
            "random_seed": 42,
            "thread_count": 12,
            "use_best_model": True,
            # "rsm":0.85,
        }

        self.model = cb.CatBoostClassifier(**params)
        self.model.fit(train_pool, eval_set=val_pool)

        val_score = self.model.predict_proba(X_val)
        y_val_flat = y_val.ravel()
        self._eval(y_val_flat, val_score, prefix="VAL")
        return self

    def _eval(self, y_true, score, prefix="TEST"):
        pred_idx = np.argmax(score, axis=1)
        print(f"\n==== {prefix} METRICS ====")
        print("\nConfusion Matrix:\n", confusion_matrix(y_true, pred_idx))
        print("\nClassification Report:\n", classification_report(y_true, pred_idx, digits=4))

    def predict_proba(self, X):
        return self.model.predict_proba(X)

    def predict_label(self, X):
        proba = self.predict_proba(X)
        return np.argmax(proba, axis=1)

    def topk(self, X, k=100, target_class=1):
        proba = self.predict_proba(X)
        class_score = proba[:, target_class]
        idx = np.argsort(class_score)[::-1][:k]
        return idx, class_score[idx]

def nn_feature_data():
    train_sample=[]
    train_labels=[]
    test_sample=[]
    test_labels=[]

    df_A = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/tgs_bench_based_random_interval/nn_random_interval_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    df_B = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/bench.v3/proof_tgs_bench/nn_tgs_bench_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    df_mix = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/bench.v3/one_bp_precise_interval/nn_one_bp_precise_interval_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)

    df_C = df_mix[df_mix['left_label'] == 1]
    df_D = df_mix[df_mix['right_label'] == 1]

    s_df_A = df_A.drop(f_list_drop1, axis=1)
    s_df_B = df_B.drop(f_list_drop1, axis=1)
    s_df_C = df_C.drop(f_list_drop2, axis=1)
    s_df_D = df_D.drop(f_list_drop2, axis=1)

    A_features = s_df_A.to_numpy()
    A_labels = np.array([0]*len(s_df_A))
    B_features = s_df_B.to_numpy()
    B_labels = np.array([1]*len(s_df_B))
    C_features = s_df_C.to_numpy()
    C_labels = np.array([2]*len(s_df_C))
    D_features = s_df_D.to_numpy()
    D_labels = np.array([3]*len(s_df_D))

    train_sample.append(A_features)
    train_labels.append(A_labels)

    df_random_large_cnv = pd.read_csv("/home/data/user/gyzhe/data/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/random.large.CNV/RP_SR_depth_BAF_proofStatistic.csv",index_col=0)
    s_df_random_large_cnv = df_random_large_cnv.drop(f_list_drop1, axis=1)
    s_df_random_large_cnv_features = s_df_random_large_cnv.to_numpy()
    s_df_random_large_cnv_labels = np.array([0] * len(s_df_random_large_cnv))
    train_sample.append(s_df_random_large_cnv_features)
    train_labels.append(s_df_random_large_cnv_labels)

    test_dataset_random = "/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/test_dataset_random_interval/nn_random_interval_RP_SR_depth_BAF_proofStatistic.csv"
    df_random = pd.read_csv(test_dataset_random,index_col=0)
    t_df_random = df_random.drop(f_list_drop1, axis=1)
    t_df_random_features = t_df_random.to_numpy()
    t_df_random_labels = np.array([0] * len(t_df_random))
    test_sample.append(t_df_random_features)
    test_labels.append(t_df_random_labels)
    test_sample.append(B_features)
    test_labels.append(B_labels)
    test_sample.append(C_features)
    test_labels.append(C_labels)
    test_sample.append(D_features)
    test_labels.append(D_labels)
    #loh_csv
    loh_csv = "/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/bed/both.bp.shift/nn_RP_SR_depth_BAF_proofStatistic.csv"
    df_loh = pd.read_csv(loh_csv,index_col=0)

    s_df_loh = df_loh.drop(f_list_drop1, axis=1)
    loh_features = s_df_loh.to_numpy()
    loh_labels = np.array([5] * len(s_df_loh))
    train_sample.append(loh_features)
    train_labels.append(loh_labels)
    loh_csv = "/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/vcf/both.bp.shift/nn_RP_SR_depth_BAF_proofStatistic.csv"
    df_loh = pd.read_csv(loh_csv,index_col=0)

    s_df_loh = df_loh.drop(f_list_drop1, axis=1)
    loh_features = s_df_loh.to_numpy()
    loh_labels = np.array([5] * len(s_df_loh))
    train_sample.append(loh_features)
    train_labels.append(loh_labels)

    loh_csv = "/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/bed/both.bp.shift10k/nn_RP_SR_depth_BAF_proofStatistic.csv"
    df_loh = pd.read_csv(loh_csv,index_col=0)
    s_df_loh = df_loh.drop(f_list_drop1, axis=1)
    loh_features = s_df_loh.to_numpy()
    loh_labels = np.array([5] * len(s_df_loh))
    test_sample.append(loh_features)
    test_labels.append(loh_labels)
    loh_csv = "/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/vcf/both.bp.shift10k/nn_RP_SR_depth_BAF_proofStatistic.csv"
    df_loh = pd.read_csv(loh_csv, index_col=0)
    s_df_loh = df_loh.drop(f_list_drop1, axis=1)
    loh_features = s_df_loh.to_numpy()
    loh_labels = np.array([5] * len(s_df_loh))
    test_sample.append(loh_features)
    test_labels.append(loh_labels)
    sim_df = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/both.bp.shift/vaf0.3.50x_nn_training_dataset_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)

    s_df = sim_df.drop(f_list_drop1, axis=1)
    sim_features = s_df.to_numpy()
    sim_labels = np.array([5] * len(s_df))
    train_sample.append(sim_features)
    train_labels.append(sim_labels)
    for vaf in ['vaf0.1.50x', 'vaf0.3.50x']:
        sim_df_B = pd.read_csv(
            '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/bench_proof/nn_csv/'+vaf+'_nn_training_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)
        sim_df_mix = pd.read_csv('/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/one_bp_precise_interval_proof/nn_csv/'+vaf+'_nn_training_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)

        sim_df_C = sim_df_mix[sim_df_mix['left_label'] == 1]
        sim_df_D = sim_df_mix[sim_df_mix['right_label'] == 1]
        s_df_C = sim_df_C.drop(f_list_drop2, axis=1)
        s_df_D = sim_df_D.drop(f_list_drop2, axis=1)

        s_df_B = sim_df_B.drop(f_list_drop1, axis=1)
        sim_B_features = s_df_B.to_numpy()
        sim_B_labels = np.array([1]*len(s_df_B))
        sim_C_features = s_df_C.to_numpy()
        sim_C_labels = np.array([2] *len(s_df_C))
        sim_D_features = s_df_D.to_numpy()
        sim_D_labels = np.array([3] *len(s_df_D))
        train_sample.append(sim_B_features)
        train_labels.append(sim_B_labels)
        train_sample.append(sim_C_features)
        train_labels.append(sim_C_labels)
        train_sample.append(sim_D_features)
        train_labels.append(sim_D_labels)

    df_HG008_bed_bench = pd.read_csv("/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/bed/proof/nn_RP_SR_depth_BAF_proofStatistic.csv",
        index_col=0)
    df_HG008_bed_mix = pd.read_csv('/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/bed/one_bp_precise_interval/nn_one_bp_precise_interval_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)

    df_C = df_HG008_bed_mix[df_HG008_bed_mix['left_label'] == 1]
    df_D = df_HG008_bed_mix[df_HG008_bed_mix['right_label'] == 1]

    s_df_B = df_HG008_bed_bench.drop(f_list_drop1, axis=1)
    s_df_C = df_C.drop(f_list_drop2, axis=1)
    s_df_D = df_D.drop(f_list_drop2, axis=1)

    B_features = s_df_B.to_numpy()
    B_labels = np.array([1] * len(s_df_B))
    C_features = s_df_C.to_numpy()
    C_labels = np.array([2] * len(s_df_C))
    D_features = s_df_D.to_numpy()
    D_labels = np.array([3] * len(s_df_D))
    target_n = 450
    repeat_factor = target_n // len(s_df_B) + 1
    B_features_oversampled = np.repeat(B_features, repeat_factor, axis=0)[:target_n]
    B_labels_oversampled = np.repeat(B_labels, repeat_factor, axis=0)[:target_n]
    C_features_oversampled = np.repeat(C_features, repeat_factor, axis=0)[:target_n]
    C_labels_oversampled = np.repeat(C_labels, repeat_factor, axis=0)[:target_n]
    D_features_oversampled = np.repeat(D_features, repeat_factor, axis=0)[:target_n]
    D_labels_oversampled = np.repeat(D_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(B_features_oversampled)
    train_labels.append(B_labels_oversampled)
    train_sample.append(C_features_oversampled)
    train_labels.append(C_labels_oversampled)
    train_sample.append(D_features_oversampled)
    train_labels.append(D_labels_oversampled)

    df_HG008_vcf_bench = pd.read_csv("/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/vcf/proof/nn_RP_SR_depth_BAF_proofStatistic.csv",
        index_col=0)
    df_HG008_vcf_mix = pd.read_csv('/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/vcf/one_bp_precise_interval/nn_one_bp_precise_interval_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    df_C = df_HG008_vcf_mix[df_HG008_vcf_mix['left_label'] == 1]
    df_D = df_HG008_vcf_mix[df_HG008_vcf_mix['right_label'] == 1]

    s_df_B = df_HG008_vcf_bench.drop(f_list_drop1, axis=1)
    s_df_C = df_C.drop(f_list_drop2, axis=1)
    s_df_D = df_D.drop(f_list_drop2, axis=1)

    B_features = s_df_B.to_numpy()
    B_labels = np.array([1] * len(s_df_B))
    C_features = s_df_C.to_numpy()
    C_labels = np.array([2] * len(s_df_C))
    D_features = s_df_D.to_numpy()
    D_labels = np.array([3] * len(s_df_D))
    target_n = 240
    repeat_factor = target_n // len(s_df_B) + 1
    B_features_oversampled = np.repeat(B_features, repeat_factor, axis=0)[:target_n]
    B_labels_oversampled = np.repeat(B_labels, repeat_factor, axis=0)[:target_n]
    C_features_oversampled = np.repeat(C_features, repeat_factor, axis=0)[:target_n]
    C_labels_oversampled = np.repeat(C_labels, repeat_factor, axis=0)[:target_n]
    D_features_oversampled = np.repeat(D_features, repeat_factor, axis=0)[:target_n]
    D_labels_oversampled = np.repeat(D_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(B_features_oversampled)
    train_labels.append(B_labels_oversampled)
    train_sample.append(C_features_oversampled)
    train_labels.append(C_labels_oversampled)
    train_sample.append(D_features_oversampled)
    train_labels.append(D_labels_oversampled)

    fake_cnv_csv = "/home/data/user/gyzhe/data/project/validateCNVBeakpoint/data/HG008/bench/fake.CNV/nn_RP_SR_depth_BAF_proofStatistic.csv"
    fake_cnv = pd.read_csv(fake_cnv_csv,index_col=0)

    neg_df = fake_cnv.drop(f_list_drop1, axis=1)

    neg_features = neg_df.to_numpy()
    neg_labels = np.array([4] * len(neg_df))
    train_sample.append(neg_features)
    train_labels.append(neg_labels)

    features = np.vstack(tuple(train_sample))
    labels = np.hstack(tuple(train_labels))
    test_feature = np.vstack(tuple(test_sample))
    test_l=np.hstack(tuple(test_labels))
    feature_names = neg_df.columns.tolist()
    return features, labels, test_feature, test_l,feature_names
def train():
    X_train_val, y_train_val, X_test, y_test, feature_names = nn_feature_data()
    model = CNVCatBoost()
    model.train(X_train_val, y_train_val)
    dir_model = "~/"
    joblib.dump(model, dir_model + "multi_catboost_model.pkl")

def model_val(model_dir,independent_CNV_callset_feature,output_csv):
    loaded_cb = joblib.load(model_dir)
    df_B = pd.read_csv(independent_CNV_callset_feature, index_col=0)
    s_df = df_B.drop(f_list_drop1, axis=1)
    test_feature = s_df.to_numpy()
    pred_label = loaded_cb.predict_label(test_feature)
    df_B['pred_class'] = pred_label
    df_B.to_csv(output_csv, index=True)


