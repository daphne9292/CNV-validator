import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import recall_score, precision_score, f1_score, roc_curve, auc
from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn.preprocessing import RobustScaler


class SimpleModel(nn.Module):
    def __init__(self, input_size):
        super(SimpleModel, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 3)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        return self.sigmoid(x)


def train(model, X_train_val, y_train_val, optimizer, num_epochs=100,val_size=0.2,best_metric='f1'):
    best_f1 = 0.0  
    best_model_state = None  #
    label_counts = y_train_val.sum(dim=0)  
    # print("label_counts: ", label_counts)
    total = y_train_val.shape[0]
    # pos_weight = (total - label_counts) / (label_counts + 1e-6) 
    # pos_weight = pos_weight.to(device)
    for epoch in range(num_epochs):
        model.train()
        train_indices = np.random.choice(len(X_train_val), size=int(len(X_train_val) * (1 - val_size)), replace=False)
        val_indices = np.setdiff1d(np.arange(len(X_train_val)), train_indices)

        X_train, y_train = X_train_val[train_indices], y_train_val[train_indices]

        X_val, y_val = X_train_val[val_indices], y_train_val[val_indices]

    
        optimizer.zero_grad()
        outputs = model(X_train).squeeze()
      
        # criterion = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        criterion = torch.nn.BCEWithLogitsLoss()
        loss = criterion(outputs, y_train)
        # loss = nn.BCELoss()(outputs, y_train)
        # loss = cross_entropy_loss(outputs, y_train)

        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val).squeeze()
            val_preds = (val_outputs > 0.5).float()
            recall = recall_score(y_val.cpu(), val_preds.cpu(),average='macro')
            precision = precision_score(y_val.cpu(), val_preds.cpu(),average='macro')
            f1 = f1_score(y_val.cpu(), val_preds.cpu(),average='macro')

        n_labels = y_val.shape[1]
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item():.4f}, Recall: {recall:.4f}, Precision: {precision:.4f}, F1: {f1:.4f}')
        for i in range(n_labels):
            fpr[i], tpr[i], _ = roc_curve(y_val[:, i].cpu(), val_outputs[:, i].cpu())
            roc_auc[i] = auc(fpr[i], tpr[i])
            print(f'Lable [{i}], AUC: {roc_auc[i]:.4f}')

        if best_metric == 'f1' and f1 > best_f1:
            best_f1 = f1
            best_model_state = model.state_dict()

    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return model  
def training_dataset_split():
    train_sample=[]
    train_labels=[]
    test_sample=[]
    label_list=[]
    # simulated CNV sample
    for vaf in['vaf0.005.100x', 'vaf0.01.100x', 'vaf0.02.100x','vaf0.1.100x', 'vaf0.3.100x', 'vaf0.5.100x', \
                'vaf0.005.50x','vaf0.01.50x', 'vaf0.02.50x', 'vaf0.1.50x', 'vaf0.3.50x', 'vaf0.5.50x', \
               'vaf0.005.20x', 'vaf0.01.20x', 'vaf0.02.20x', 'vaf0.1.20x', 'vaf0.3.20x', 'vaf0.5.20x',\
                'vaf0.005.10x','vaf0.01.10x', 'vaf0.02.10x', 'vaf0.1.10x', 'vaf0.3.10x', 'vaf0.5.10x',\
               'vaf0.005.5x','vaf0.01.5x', 'vaf0.02.5x', 'vaf0.1.5x', 'vaf0.3.5x', 'vaf0.5.5x']:
        # positive sample [1,1,1]
        sim_df_B = pd.read_csv(
            '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/bench_proof/'+vaf+'_training_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)
        # negative sample [0,0,1] [0,1,0]
        sim_df_mix = pd.read_csv(
            '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/one_bp_precise_interval_proof/'+vaf+'_training_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)
        sim_df_C = sim_df_mix[sim_df_mix['left_label'] == 1]
        sim_df_D = sim_df_mix[sim_df_mix['right_label'] == 1]

        sim_df_B.fillna(0, inplace=True)  
        sim_df_C.fillna(0, inplace=True)  
        sim_df_D.fillna(0, inplace=True) 
        # select some feature for training
        s_df_B = sim_df_B.drop(['svtype', 'report_len','valid_t_SR_percent_left_bp_interval','valid_t_SR_percent_right_bp_interval','valid_n_SR_percent_left_bp_interval','valid_n_SR_percent_right_bp_interval'], axis=1)
        s_df_C = sim_df_C.drop(['svtype', 'left_label', 'right_label', 'report_len','valid_t_SR_percent_left_bp_interval','valid_t_SR_percent_right_bp_interval','valid_n_SR_percent_left_bp_interval','valid_n_SR_percent_right_bp_interval'], axis=1)
        s_df_D = sim_df_D.drop(['svtype', 'left_label', 'right_label', 'report_len','valid_t_SR_percent_left_bp_interval','valid_t_SR_percent_right_bp_interval','valid_n_SR_percent_left_bp_interval','valid_n_SR_percent_right_bp_interval'], axis=1)

        sim_B_features = s_df_B.to_numpy()
        sim_B_labels = np.array([[1, 1, 1] for _ in range(s_df_B.shape[0])])
        sim_C_features = s_df_C.to_numpy()
        sim_C_labels = np.array([[0, 1, 0] for _ in range(s_df_C.shape[0])])
        sim_D_features = s_df_D.to_numpy()
        sim_D_labels = np.array([[0, 0, 1] for _ in range(s_df_D.shape[0])])

        train_sample.append(sim_B_features)
        train_labels.append(sim_B_labels)
        train_sample.append(sim_C_features)
        train_labels.append(sim_C_labels)
        train_sample.append(sim_D_features)
        train_labels.append(sim_D_labels)

    # negative sample[0,0,0],random selected region in Chinese Quartet genome
    sim_df_A = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/random_interval_proof/RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    sim_df_A.fillna(0, inplace=True)  
    # select some feature for training
    s_df_A = sim_df_A.drop(['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    sim_A_features = s_df_A.to_numpy()
    sim_A_labels = np.array([[0, 0, 0] for _ in range(s_df_A.shape[0])])
    # upsample negative sample
    target_n = 3000
    repeat_factor = target_n // len(sim_A_features) + 1
    A_features_oversampled = np.repeat(sim_A_features, repeat_factor, axis=0)[:target_n]
    A_labels_oversampled = np.repeat(sim_A_labels, repeat_factor, axis=0)[:target_n]

    train_sample.append(A_features_oversampled)
    train_labels.append(A_labels_oversampled)
    # negative sample [0,0,0],cnvkit call on ChineseQuartet LCL5(as normal)/LCL6(as tumor) 2300
    cnvkit_negative = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/cnvkit/LCL5normal/proof/RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    cnvkit_negative.fillna(0, inplace=True)
    neg_df = cnvkit_negative.drop(
        ['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    neg_features = neg_df.to_numpy()
    neg_labels = np.array([[0, 0, 0] for _ in range(neg_df.shape[0])])
    target_n = 4000
    repeat_factor = target_n // len(neg_features) + 1
    neg_features_oversampled = np.repeat(neg_features, repeat_factor, axis=0)[:target_n]
    neg_labels_oversampled = np.repeat(neg_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(neg_features_oversampled)
    train_labels.append(neg_labels_oversampled)
    # negative sample [0,0,0],cnvkit call on ChineseQuartet LCL6(as normal)/LCL5(as tumor)
    cnvkit_negative = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/cnvkit/LCL6normal/proof/RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    cnvkit_negative.fillna(0, inplace=True)  # 用0填充缺失值
    neg_df = cnvkit_negative.drop(['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    neg_features = neg_df.to_numpy()
    neg_labels = np.array([[0, 0, 0] for _ in range(neg_df.shape[0])])
    target_n = 5000
    repeat_factor = target_n // len(neg_features) + 1
    neg_features_oversampled = np.repeat(neg_features, repeat_factor, axis=0)[:target_n]
    neg_labels_oversampled = np.repeat(neg_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(neg_features_oversampled)
    train_labels.append(neg_labels_oversampled)
    # negative sample [0,0,0],facets call on ChineseQuartet LCL5(as normal)/LCL6(as tumor)
    facets_negative = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/facets/LCL5normal/proof/RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    facets_negative.fillna(0, inplace=True) 
    neg_df = facets_negative.drop(['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    neg_features = neg_df.to_numpy()
    neg_labels = np.array([[0, 0, 0] for _ in range(neg_df.shape[0])])
    target_n = 2500
    repeat_factor = target_n // len(neg_features) + 1
    neg_features_oversampled = np.repeat(neg_features, repeat_factor, axis=0)[:target_n]
    neg_labels_oversampled = np.repeat(neg_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(neg_features_oversampled)
    train_labels.append(neg_labels_oversampled)
    # negative sample [0,0,0],sequenza call on ChineseQuartet LCL5(as normal)/LCL6(as tumor)
    sequenza_negative = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/data/ChineseQuartet/analysis/methods.call.cnv/sequenza/LCL5normal/proof/RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    sequenza_negative.fillna(0, inplace=True) 
    neg_df = sequenza_negative.drop(['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    neg_features = neg_df.to_numpy()
    neg_labels = np.array([[0, 0, 0] for _ in range(neg_df.shape[0])])
    target_n = 2500
    repeat_factor = target_n // len(neg_features) + 1
    neg_features_oversampled = np.repeat(neg_features, repeat_factor, axis=0)[:target_n]
    neg_labels_oversampled = np.repeat(neg_labels, repeat_factor, axis=0)[:target_n]
    train_sample.append(neg_features_oversampled)
    train_labels.append(neg_labels_oversampled)
    features = np.vstack(tuple(train_sample))
    labels = np.vstack(tuple(train_labels))

    # test set
    # positive sample,HCC1395 benchmark,label[1,1,1]
    df_B = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/proof_tgs_bench/tgs_bench_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)
    df_B.fillna(0, inplace=True) 
    # negative sample,shifted breakpoint of HCC1395 benchmark,label[0,0,1] [0,1,0]
    df_mix = pd.read_csv(
        '/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/one_bp_precise_interval/one_bp_precise_interval_RP_SR_depth_BAF_proofStatistic.csv',
        index_col=0)

    df_C = df_mix[df_mix['left_label'] == 1]
    df_D = df_mix[df_mix['right_label'] == 1]
    df_C.fillna(0, inplace=True)  
    df_D.fillna(0, inplace=True)  

    s_df_B = df_B.drop(['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
         'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
    s_df_C = df_C.drop(['svtype', 'left_label', 'right_label', 'report_len', 'valid_t_SR_percent_left_bp_interval',
                        'valid_t_SR_percent_right_bp_interval', 'valid_n_SR_percent_left_bp_interval',
                        'valid_n_SR_percent_right_bp_interval'], axis=1)
    s_df_D = df_D.drop(['svtype', 'left_label', 'right_label', 'report_len', 'valid_t_SR_percent_left_bp_interval',
                        'valid_t_SR_percent_right_bp_interval', 'valid_n_SR_percent_left_bp_interval',
                        'valid_n_SR_percent_right_bp_interval'], axis=1)
    B_features = s_df_B.to_numpy()
    B_labels = np.array([[1, 1, 1] for _ in range(s_df_B.shape[0])])
    C_features = s_df_C.to_numpy()
    C_labels = np.array([[0, 1, 0] for _ in range(s_df_C.shape[0])])
    D_features = s_df_D.to_numpy()
    D_labels = np.array([[0, 0, 1] for _ in range(s_df_D.shape[0])])
    test_sample.append(B_features)
    test_sample.append(C_features)
    test_sample.append(D_features)
    label_list.append(B_labels)
    label_list.append(C_labels)
    label_list.append(D_labels)
    for vaf in ['vaf0.005.100x', 'vaf0.01.100x', 'vaf0.02.100x', 'vaf0.1.100x', 'vaf0.3.100x', 'vaf0.5.100x', \
                'vaf0.005.50x', 'vaf0.01.50x', 'vaf0.02.50x', 'vaf0.1.50x', 'vaf0.3.50x', 'vaf0.5.50x', \
                'vaf0.005.20x', 'vaf0.01.20x', 'vaf0.02.20x', 'vaf0.1.20x', 'vaf0.3.20x', 'vaf0.5.20x', \
                'vaf0.005.10x', 'vaf0.01.10x', 'vaf0.02.10x', 'vaf0.1.10x', 'vaf0.3.10x', 'vaf0.5.10x', \
                'vaf0.005.5x', 'vaf0.01.5x', 'vaf0.02.5x', 'vaf0.1.5x', 'vaf0.3.5x', 'vaf0.5.5x']:
        test_B = pd.read_csv(
            '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/bench_proof/' + vaf + '_test_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)
        test_mix = pd.read_csv(
            '/data/home/gyzhe/project/validateCNVBeakpoint/sim/two.abnormal.size/one_bp_precise_interval_proof/' + vaf + '_test_dataset_RP_SR_depth_BAF_proofStatistic.csv',
            index_col=0)
        test_C = test_mix[test_mix['left_label'] == 1]
        test_D = test_mix[test_mix['right_label'] == 1]
        test_B.fillna(0, inplace=True)  
        test_C.fillna(0, inplace=True)  
        test_D.fillna(0, inplace=True)  

        t_df_B = test_B.drop(
            ['svtype', 'report_len', 'valid_t_SR_percent_left_bp_interval', 'valid_t_SR_percent_right_bp_interval',
             'valid_n_SR_percent_left_bp_interval', 'valid_n_SR_percent_right_bp_interval'], axis=1)
        t_df_C = test_C.drop(
            ['svtype', 'left_label', 'right_label', 'report_len', 'valid_t_SR_percent_left_bp_interval',
             'valid_t_SR_percent_right_bp_interval', 'valid_n_SR_percent_left_bp_interval',
             'valid_n_SR_percent_right_bp_interval'], axis=1)
        t_df_D = test_D.drop(
            ['svtype', 'left_label', 'right_label', 'report_len', 'valid_t_SR_percent_left_bp_interval',
             'valid_t_SR_percent_right_bp_interval', 'valid_n_SR_percent_left_bp_interval',
             'valid_n_SR_percent_right_bp_interval'], axis=1)

        t_B_features = t_df_B.to_numpy()
        t_B_labels = np.array([[1, 1, 1] for _ in range(t_df_B.shape[0])])
        t_C_features = t_df_C.to_numpy()
        t_C_labels = np.array([[0, 1, 0] for _ in range(t_df_C.shape[0])])
        t_D_features = t_df_D.to_numpy()
        t_D_labels = np.array([[0, 0, 1] for _ in range(t_df_D.shape[0])])
        test_sample.append(t_B_features)
        test_sample.append(t_C_features)
        test_sample.append(t_D_features)
        label_list.append(t_B_labels)
        label_list.append(t_C_labels)
        label_list.append(t_D_labels)
    test_features = np.vstack(test_sample)
    test_labels = np.vstack(label_list)
    return features, labels, test_features, test_labels

def model_train(device):
    X_train_val,y_train_val,X_test,y_test = training_dataset_split()

    # features standardized
    manual_scaler = RobustScaler()
    manual_scaler.fit(X_train_val)
    features_scaled = manual_scaler.transform(X_train_val)
    test_scaled = manual_scaler.transform(X_test)

    # turn to tensor
    X_train_val = torch.tensor(features_scaled, dtype=torch.float32).to(device)
    y_train_val = torch.tensor(y_train_val, dtype=torch.float32).to(device)
    X_test = torch.tensor(test_scaled, dtype=torch.float32).to(device)

    # model parameters  initialized
    model = SimpleModel(X_train_val.shape[1]).to(device)
    # optimizer
    optimizer = optim.Adam(model.parameters(), lr=0.002)

    best_model = train(model, X_train_val, y_train_val, optimizer, num_epochs=1200)
    # save model
    torch.save(best_model.state_dict(), "/data/home/gyzhe/project/validateCNVBeakpoint/output.hcc1395/two.abnormal.size/supervised_learning/8/trained_model.pth")

    model = best_model
    # model evaluation
    model.eval()
    with torch.no_grad():
        test_outputs = model(X_test).squeeze()
        test_preds = (test_outputs > 0.5).float()
        recall = recall_score(y_test, test_preds.cpu(), average='macro')
        precision = precision_score(y_test, test_preds.cpu(), average='macro')
        f1 = f1_score(y_test, test_preds.cpu(), average='macro')
        n_labels = y_test.shape[1]
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        all_aucs = []
        for i in range(n_labels):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], test_outputs[:, i].cpu())
            roc_auc[i] = auc(fpr[i], tpr[i])
            all_aucs.append(roc_auc[i])
            print(f'Test Set Lable [{i}], AUC: {roc_auc[i]:.4f}')

        # macro AUC
        macro_auc = sum(all_aucs) / len(all_aucs)
        print(f'Macro Average AUC: {macro_auc:.4f}')
    print(f'Test Set Recall: {recall:.4f}')
    print(f'Test Set Precision: {precision:.4f}')
    print(f'Test Set F1: {f1:.4f}')

    
    # plt.figure()
    # for i in range(n_labels):
    #     plt.plot(fpr[i], tpr[i], label='ROC curve of label {0} (area = {1:0.2f})'
    #                                    ''.format(i, roc_auc[i]))
    #
    # plt.plot([0, 1], [0, 1], 'k--')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('ROC Curve for Each Label')
    # plt.legend(loc="lower right")
    # plt.savefig('./test_set_roc.png')
def model_val(device,model_dir,independent_CNV_callset_feature,output_csv):
    X_train_val,y_train_val,X_test,y_test = training_dataset_split()

    manual_scaler = RobustScaler()
    manual_scaler.fit(X_train_val)

    loaded_model = SimpleModel(X_train_val.shape[1])
    loaded_model.load_state_dict(torch.load(model_dir))
    loaded_model.to(device)

    loaded_model.eval()

    df = pd.read_csv(independent_CNV_callset_feature,index_col=0)
    df.fillna(0, inplace=True)  
    s_df = df.drop(['svtype', 'report_len','valid_t_SR_percent_left_bp_interval','valid_t_SR_percent_right_bp_interval','valid_n_SR_percent_left_bp_interval','valid_n_SR_percent_right_bp_interval'], axis=1)
    test_feature = s_df.to_numpy()
    test_scaled = manual_scaler.transform(test_feature)
    X_test = torch.tensor(test_scaled, dtype=torch.float32).to(device)
    with torch.no_grad():
        test_outputs = loaded_model(X_test).squeeze()
        test_preds = (test_outputs > 0.5).float()
    df['pred_cnv'] = test_preds.cpu()[:, 0]
    df['pred_lbp_label'] = test_preds.cpu()[:, 1]
    df['pred_rbp_label'] = test_preds.cpu()[:, 2]
    df.to_csv(output_csv, index=True)





