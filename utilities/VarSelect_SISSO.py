# Created by Zhen Guo and Runhai Ouyang, 2022.2
# Variable Selection for SISSO (VS-SISSO), J. Chem. Theory Comput. 18, 4945 (2022)
# Usage: prepare the SISSO.in and train.dat in working directory, and then run the program
# by using python3 or newer versions,e.g.:  python3 VarSelect_SISSO.py

# version_20250320 fix a bug : the number of features in train.dat is not the same as the number of features in SISSO.in

############################################ User Input #####################################
# Set your control parameters here
n_init = 10           # Initial size of the S in [J. Chem. Theory Comput. 18, 4945 (2022)]
n_RS = 4              # The size of the Sa 
n_max =20             # Maximal size of the S (S=Sa+Sb)
nstep_max =100       # Maximal iterations                                                      
nstep_converge = 20   # Converged if the best model unchanged after given number of steps.   
restart = 0           # 0: start from scratch, 1: continue the unfinished job
runSISSO = 'mpirun -np 64  SISSO.3.5 > SISSO.log'  # The command to run SISSO
#############################################################################################




import os
import copy
import time
import random
import math


def SISSO_out_reader(SISSO_out_folder, dimension, all_features, maths_operators):
    # From SISSO.out, read the errors and involved primary features in all the models

    SISSO_out_file = open('%s/SISSO.out' % SISSO_out_folder, 'r').readlines()
    for i in range(len(SISSO_out_file)):
        SISSO_out_file[i] = SISSO_out_file[i].strip()
    feature_list = []
    descriptor_dict = {}
#    descriptor_dim = 0
    descriptor_dim = dimension
    for i in range(len(SISSO_out_file)):
        if SISSO_out_file[i].count('D descriptor'):
#            descriptor_dim += 1
            descriptor_list = []

            for d in range(descriptor_dim):
                descriptor = SISSO_out_file[i + d + 2]

                for st in range(len(descriptor)):
                    if descriptor[st] == '=':
                        descriptor = descriptor[st + 1:].split()[0]
                        break
                for feature in descriptor_2_features(descriptor, all_features, maths_operators):
                    feature_list.append(feature)
                feature_list = list(set(feature_list))
                descriptor_list.append(descriptor)
            descriptor_dict.update({descriptor_dim: feature_list})
    total_RMSE_dict = {}
#    RMSE_dim = 0
    RMSE_dim = dimension
    for i in range(len(SISSO_out_file)):
        if SISSO_out_file[i].count('RMSE and MaxAE:'):
#            RMSE_dim += 1
            RMSE_str = SISSO_out_file[i].replace('RMSE and MaxAE:', '').strip()
            RMSE = ''
            for j in list(RMSE_str):
                if j != ' ':
                    RMSE += j
                else:
                    break
            total_RMSE_dict.update({RMSE_dim: float(RMSE)})
    return feature_list, total_RMSE_dict[dimension]

def descriptor_2_features(descriptor, all_features,maths_operators):
# Identify the primary features in a descriptor formula
    import copy

    brace = []
    brace_position = []
    for i in range(len(descriptor)):
        if descriptor[i] == '(':
            brace.append(0)
            brace_position.append(i)
        if descriptor[i] == ")":
            brace.append(1)
            brace_position.append(i)

    features = []

    while brace:
        for i in range(len(brace)):
            if (brace[i] == 0) and (brace[i + 1] == 1):
                features.append(descriptor[brace_position[i] + 1:brace_position[i + 1]])
                # if features[-1].startswith('('):
                #     del features[-1]

                del brace[i:i + 2]
                del brace_position[i:i + 2]
                break

    features_new = []
    for feature in features:
        features_new.append(feature)

    for Feature in features:
        maths_operator_position = []
        maths_operator_length = []
        for i in range(len(Feature)):
            for operator in maths_operators:
                op_len = len(operator)
                if Feature[i:i + op_len] == operator:
                    maths_operator_position.append(i)
                    maths_operator_length.append(op_len)
                    break
        Feature_cp = copy.copy(Feature)
        count = 0
        count_max = len(copy.copy(maths_operator_position))
        while count < count_max:

            for j in range(len(maths_operator_position)):
                features_new.append(Feature_cp[:maths_operator_position[j]])
                features_new.append(Feature_cp[maths_operator_position[j] + maths_operator_length[j]:])

            maths_operator_length_0 = maths_operator_length[:1][0] + maths_operator_position[:1][0]
            Feature_cp = Feature_cp[maths_operator_length_0:]
            del maths_operator_length[:1]
            del maths_operator_position[:1]
            for j in range(len(maths_operator_position)):
                maths_operator_position[j] = maths_operator_position[j] - maths_operator_length_0

            count += 1

    features_out = []
    for i in features_new:
        if (i not in features_out) & (i in all_features):
            features_out.append(i)

    return features_out

def initial_SISSO_in_2_output_parameter(initial_file_dir, all_features_list, output_parameter):
    # Read information from SISSO.in

    parameter_startwith_dict = {1: "funit", 2: "ops", 3: "desc_dim", 4:"nsf"}

    SISSO_in_file = open('%s/SISSO.in' % initial_file_dir, 'r').readlines()
    for i in range(len(SISSO_in_file)):
        SISSO_in_file[i] = SISSO_in_file[i].strip()
    output_para_in_file = ''
    for i in range(len(SISSO_in_file)):
        if SISSO_in_file[i].startswith(parameter_startwith_dict[output_parameter]):
            output_para_in_file = SISSO_in_file[i]
            # output_para_in_file = output_para_in_file.replace(parameter_startwith_dict[output_parameter], '')
            output_para_in_file = output_para_in_file.split('=')[1]
    if output_parameter == 1:  # funit
        funit = output_para_in_file
        funit = funit.replace('(', "")
        for alp in range(len(funit)):
            if funit[alp] == '!':
                funit = funit[:alp].strip()
                break
        if funit == ')':
            return [[]]
        else:
            funit_list = []
            operator = ''
            for i in funit:
                if i == '!':
                    break
                if (i != ':' and i != ")"):
                    operator += i
                else:
                    funit_list.append(int(operator))
                    operator = ''

            if (funit_list[-1] > len(all_features_list)):
                exit('Error: funit out of range！\n'
                     'Check the parameter \'funit\' and \'nsf\' in SISSO.in')
            if len(funit_list) % 2 != 0:
                exit('Error: wrong \'funit\' setting! \n'
                     'Check the parameter \'funit\' in SISSO.in')
            feature_class = []
            list_new = []
            for i in range(len(funit_list)):
                if i % 2 == 1:
                    if i != len(funit_list) - 1:
                        list_new += all_features_list[funit_list[i]:funit_list[i + 1] - 1]
                    else:
                        list_new += all_features_list[funit_list[i]:]
                    continue
                feature_class.append(all_features_list[funit_list[i] - 1:funit_list[i + 1]])
            # if list_new != []:
            #     feature_class.append(list_new)
            return feature_class


    elif output_parameter == 2:  # ops
        operators_list = []
        operators = output_para_in_file.split('!')[0].strip()
        operators = operators.replace("','", "")
        operators = operators.replace(")(", " ")
        operators = operators.replace("'(", "").replace(")'", "")
        operators_list = list(set(operators.split()))
        return operators_list
    elif output_parameter == 3:
        desc_dim = output_para_in_file
        desc_dim = desc_dim.replace('desc_dim=', '')
        desc_dim = desc_dim[:3]
        desc_dim = desc_dim.strip()
        return int(desc_dim)
    elif output_parameter == 4: 
        nsf = output_para_in_file
        nsf = nsf.split("nsf=")[-1].split("!")[0].strip()
        return int(nsf)

def initial_train_dat_2_output_parameter(train_dat_folder, output_parameter):
    # Read data from train.dat

    train_dat_lines = open('%s/train.dat' % train_dat_folder).readlines()
    for line in range(len(train_dat_lines)):
        train_dat_lines[line] = train_dat_lines[line].replace(',',' ').replace('\t',' ')
        train_dat_lines[line] = train_dat_lines[line].strip()
        if not train_dat_lines[line]:
            train_dat_lines.remove('')
    for line in range(len(train_dat_lines)):
        train_dat_lines[line] = train_dat_lines[line].split()
    features_name_list = train_dat_lines[0][2:]
    materials_name = train_dat_lines[0][0]
    property_name = train_dat_lines[0][1]
    train_dat = {}
    for line in range(len(train_dat_lines)):
        if line == 0:
            for num in range(len(train_dat_lines[line])):
                train_dat.update({train_dat_lines[line][num]: []})
        else:
            for num in range(len(train_dat_lines[line])):
                list_temp = train_dat[train_dat_lines[0][num]]
                list_temp.append(train_dat_lines[line][num])
                train_dat.update({train_dat_lines[0][num]: list_temp})

    if output_parameter == 1:
        return materials_name
    elif output_parameter == 2:
        return property_name
    elif output_parameter == 3:
        return list(features_name_list)
    elif output_parameter == 4:
        return train_dat


def build_SISSO_in(initial_SISSO_in_folder, new_SISSO_in_folder, new_features_class, features_list):
    # Update SISSO.in for new iteration
    import os
    # new_features_class = [ ["f1","f2","f3"], ["f4","f5"],["f6"] ]
    number_feature = len(features_list)
    if new_features_class == []:
        f_unit = 'funit=()\n'
    else:
        n_group = []
        for i in range(len(new_features_class)):
            n_group.append(len(new_features_class[i]))

        f_unit_list = [1]
        f_unit = 'funit=(1:'
        for i in range(len(n_group)):
            if i == 0:
                f_unit_list.append(f_unit_list[0] + n_group[i] - 1)
                f_unit += (str(f_unit_list[-1]) + ')')
            else:
                f_unit_list.append(f_unit_list[-1] + 1)
                f_unit += ('(' + str(f_unit_list[-1]) + ':')
                f_unit_list.append(n_group[i] - 1 + f_unit_list[-1])
                f_unit += (str(f_unit_list[-1]) + ')')
        f_unit += '\n'
    nsf = 'nsf=%s\n' % number_feature
    SISSO_in = open('%s/SISSO.in' % initial_SISSO_in_folder, 'r').readlines()
    for i in range(len(SISSO_in)):
        # SISSO_in[i] = SISSO_in[i].lstrip()
        if SISSO_in[i].startswith('funit'):
            SISSO_in[i] = f_unit
        if SISSO_in[i].startswith('nsf'):
            SISSO_in[i] = nsf

    if os.path.exists(new_SISSO_in_folder):
        open('%s/SISSO.in' % new_SISSO_in_folder, 'w').writelines(SISSO_in)
    else:
        os.mkdir(new_SISSO_in_folder)
        open('%s/SISSO.in' % new_SISSO_in_folder, 'w').writelines(SISSO_in)


def features_classification(features_list, all_features_class):
    # Group the primary features for creating new train.dat according to their dimensions/units
    features_class = []
    for i in all_features_class:
        list_new = list(set(i).intersection(features_list))
        features_class.append(list_new)

    features_class_new = []
    for i in features_class:
        if i:
            features_class_new.append(i)

    return features_class_new


def build_train_dat(new_train_dat_folder, new_features_class, initial_train_dat, compounds_column_name,
                    property_column_name, features_list):
    # Creat train.dat for new iterations.
    import copy
    dimensionless_features = copy.copy(features_list)
    if new_features_class:
        for i in new_features_class:
            for j in i:
                dimensionless_features.remove(j)
    new_train_dat_lines = []
    sample_num = len(initial_train_dat[property_column_name])
    for tmp_0 in range(sample_num):
        if tmp_0 == 0:
            tmp_line = ''
            for tmp in (compounds_column_name, property_column_name):
                tmp_line += '%s  ' % tmp
            for tmp_1 in new_features_class + [dimensionless_features]:
                for tmp_2 in tmp_1:
                    if tmp_1:
                        tmp_line += '%s  ' % tmp_2
            new_train_dat_lines.append(tmp_line + '\n')
        tmp_line = ''
        for tmp in (initial_train_dat[compounds_column_name][tmp_0], initial_train_dat[property_column_name][tmp_0]):
            tmp_line += '%s  ' % tmp
        for tmp_1 in new_features_class + [dimensionless_features]:
            for tmp_2 in tmp_1:
                if tmp_2:
                    tmp_line += '%s  ' % initial_train_dat[tmp_2][tmp_0]
        new_train_dat_lines.append(tmp_line + '\n')
    open('%s/train.dat' % new_train_dat_folder, 'w').writelines(new_train_dat_lines)


def check_done(task_folder):
    # Check if the SISSO job was done successfully.

    import os, time
    file_list = []
    for root, folders, files in os.walk(task_folder):
        for j in files:
            file_list.append(j)

    if 'SISSO.out' in file_list:
        SISSO_out_read = open('%s/SISSO.out' % task_folder, 'r').readlines()
        if len(SISSO_out_read) != 0:
            if SISSO_out_read[-4].startswith('Total time (second):'):
                os.system('rm -rf %s/feature_space' % task_folder)
                return 1
            else:
                return 0
        else:
            return 0
    else:
        return 0


def random_features_list(all_features, selected_features, alpha_dict,  n_init):
    # Update of the primary features for new train.dat

    unselected_features_list = list(set(all_features) - set(selected_features))
    rand_list = []
    if 1 in alpha_dict.values():
        for i in unselected_features_list:
            rand_list.append([i, random.random() * (alpha_dict[i])])
    if 1 not in alpha_dict.values():
        for i in unselected_features_list:
            rand_list.append([i, random.random()])
    # bubble sort
    for i in range(len(rand_list)):
        for j in range(len(rand_list) - i - 1):
            if rand_list[j][1] < rand_list[j + 1][1]:
                rand_list[j], rand_list[j + 1] = rand_list[j + 1], rand_list[j]
    feature_new = []
    for i in rand_list[: n_init - len(selected_features)]:
        feature_new.append(i[0])

    return selected_features + feature_new


def update_alpha_list(alpha_dict, selected_features, features_list, alpha):
    # Update the penalty factor

    # features_list_dropped = list(set(features_list) - set(selected_features))
    for i in features_list:
        alpha_old = alpha_dict[i]
        alpha_dict.update({i: alpha_old * alpha})

    return alpha_dict


def read_feature_list_from_train_data(task_folder):
    train_dat = open('%s/train.dat' % task_folder, 'r').readlines()
    columns = train_dat[0].split()
    return columns[2:], len(columns[2:])


def check_last_step():
    os.system('ls -F | grep \'/$\' > .temp_file')
    time.sleep(5)
    dir_list = open('./.temp_file', 'r').readlines()
    for i in range(len(dir_list)):
        dir_list[i] = dir_list[i].strip()[:-1]
    max_num = -1
    for i in dir_list:
        if i.isnumeric():
            if max_num < int(i):
                max_num = int(i)

    IF_continue_step_done = 0
    if max_num != -1:
        file_list = os.listdir('./iter%s' % str(max_num))
        if 'SISSO.out' in file_list:
            SISSO_out_read = open('%s/SISSO.out' % str(max_num), 'r').readlines()
            if len(SISSO_out_read) != 0:
                if SISSO_out_read[-4].startswith('Total time (second):'):
                    IF_continue_step_done = 1
    os.system('rm -f .temp_file')
    return max_num, IF_continue_step_done


# -----------------------------
time_start = time.time()
initial_file_folder = './'
compounds_column_name = initial_train_dat_2_output_parameter(initial_file_folder, 1)
property_column_name = initial_train_dat_2_output_parameter(initial_file_folder, 2)
all_features = initial_train_dat_2_output_parameter(initial_file_folder, 3)
train_dat = initial_train_dat_2_output_parameter(initial_file_folder, 4)

initial_maths_operators = initial_SISSO_in_2_output_parameter(initial_file_folder, all_features, 2)
num_of_features = initial_SISSO_in_2_output_parameter(initial_file_folder, all_features, 4)
all_features = all_features[:num_of_features]
all_features_class = initial_SISSO_in_2_output_parameter(initial_file_folder, all_features, 1)
desc_dim = initial_SISSO_in_2_output_parameter(initial_file_folder, all_features, 3)
selected_features = []
selected_features_list = []

features_list = []
# train_features=[]
features_list_list = []
RMSE_list = []
min_RMSE_list = []
VS_results = open('./VS_results', 'a')
VS_results.write(
    "Iter\t%visited\tRMSE_this\tRMSE_min\t[Var_preserved]\t[Var_discarded]\t\n")
VS_results.close()
alpha_dict = {}
for i in all_features:
    alpha_dict.update({i: 1})
min_RMSE_step = 0
visited_set = set()
alpha = 0

continue_step, IF_continue_step_done = check_last_step()
if continue_step >= 0:
  #  os.system('cp VS_results VS_results_old')
    if IF_continue_step_done == 0:
        os.system('rm -rf %s' % str(continue_step))

for i in range(nstep_max):
    if restart :

        if i < continue_step + IF_continue_step_done:
            features_list,  n_init = read_feature_list_from_train_data('./iter%s' % str(i))
            selected_features, RMSE = SISSO_out_reader('./iter%s' % str(i), desc_dim, all_features, initial_maths_operators)
            alpha_dict = update_alpha_list(alpha_dict, selected_features, features_list, alpha)
        # train_features = initial_train_dat_2_output_parameter('./iter%s' % str(i), 3)
            features_list_list.append(features_list)
            RMSE_list.append(RMSE)
            selected_features_list.append(selected_features)
            visited_set.update(features_list)
            completed_percent = float(len(visited_set)) / float(len(all_features))
            if RMSE < RMSE_list[min_RMSE_step]:
                min_RMSE_step = i
            elif RMSE == RMSE_list[min_RMSE_step] and len(selected_features) <= len(selected_features_list[min_RMSE_step]):
                min_RMSE_step = i
            else:
                selected_features, RMSE = SISSO_out_reader('./iter%s' % str(min_RMSE_step), desc_dim, all_features,
                                                       initial_maths_operators)
 
            min_RMSE_list.append(RMSE_list[min_RMSE_step])
            continue

    VS_results = open('./VS_results', 'a')
    VS_log = open("./VS_log", 'a')
    new_folder = './iter%s' % str(i)
    try:
        os.mkdir(new_folder)
    except:
        print('folder %s/ already exist!\n' % str(i))

    VS_log.write('==========' * 10)
    VS_log.write('\nIteration\t%s\n' % str(i))
    alpha_dict = update_alpha_list(alpha_dict, selected_features, features_list, alpha)
    features_list = random_features_list(all_features, selected_features, alpha_dict,  n_init)
    new_features_class = features_classification(features_list, all_features_class)
    build_SISSO_in(initial_file_folder, new_folder, new_features_class, features_list)
    build_train_dat(new_folder, new_features_class, train_dat, compounds_column_name, property_column_name,features_list)
    os.chdir(new_folder)
    os.system('%s' % runSISSO)
    os.chdir('../')

    # time.sleep(5)
    # while True:
    #     check_num = check_done(new_folder)
    #     time.sleep(5)
    #     if check_num == 1:
    #         break

    selected_features, RMSE = SISSO_out_reader('./iter%s' % str(i), desc_dim, all_features, initial_maths_operators)
    features_list_list.append(features_list)
    RMSE_list.append(RMSE)
    selected_features_list.append(selected_features)
    features_list = initial_train_dat_2_output_parameter('./iter%s' % str(i), 3)
    visited_set.update(features_list)
    completed_percent = float(len(visited_set)) / float(len(all_features))
    # if len(visited_set) == len(all_features):
    #     alpha = 1
    #     for f in all_features:
    #         alpha_dict.update({f: 1})

    if RMSE < RMSE_list[min_RMSE_step]:
        min_RMSE_step = i
    elif RMSE == RMSE_list[min_RMSE_step] and len(selected_features) <= len(selected_features_list[min_RMSE_step]):
        min_RMSE_step = i
    else:
        selected_features, RMSE = SISSO_out_reader('./iter%s' % str(min_RMSE_step), desc_dim, all_features,
                                                   initial_maths_operators)
    VS_results.write(' %s\t%s\t%.06f\t%s\t%s\t%s\t\n' % (
        str(i), format(completed_percent, '.1%'), RMSE_list[-1], RMSE_list[min_RMSE_step],
        selected_features, list(set(features_list).difference(set(selected_features)))))

    n_init = len(selected_features) + n_RS
    if n_init > n_max:
        VS_results.write('Warning: The subset size hits maximum, the n_RS for the next step is reduced to %s\n' % str(
            n_RS -(n_init - n_max)))
        n_init = n_max

    VS_log.write("Unvisited variables (labeled as 1) : \n%s\n\n" % alpha_dict)
    #VS_log.write('Number of variables in the lowest-RMSE-model is %s, at iteration %s, with the RMSE %06f\n' % (
    #    str(len(selected_features)), str(min_RMSE_step), RMSE_list[min_RMSE_step]))
    #VS_log.write('Number of variables for the next SISSO is %s \n' % str( n_init))
    VS_log.close()
    min_RMSE_list.append(RMSE_list[min_RMSE_step])
    if len(RMSE_list) >= nstep_converge:
        if len(list(set(min_RMSE_list[-nstep_converge:]))) == 1:
            VS_results.write('Converged! \n')
            VS_results.close()
            break
time_end = time.time()
VS_results = open('./VS_results', 'a')
VS_results.write('%s\t%s\n' % ('Total time (second):', str(round(time_end - time_start, 2))))
VS_results.close()
