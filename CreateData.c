#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main() {
    // --- 数据池定义 ---
    const char *doctor_names[] = {"Dr. Smith", "Dr. Jones", "Dr. Li", "Dr. Garcia", "Dr. Wang", "Dr. Chen"};
    const char *genders[] = {"Male", "Female"};
    const char *departments[] = {"Cardiology", "Neurology", "Oncology", "Pediatrics", "Orthopedics", "General"};
    const char *diagnoses[] = {"Hypertension", "Common Cold", "Influenza", "Migraine", "Asthma", "Diabetes Type 2", "Sprain"};
    const char *treatments[] = {"Medication", "Physical Therapy", "Rest and Fluids", "Surgical Procedure", "Lifestyle Changes", "Observation"};

    // 计算每个数据池的大小
    int num_doctors = sizeof(doctor_names) / sizeof(doctor_names[0]);
    int num_genders = sizeof(genders) / sizeof(genders[0]);
    int num_departments = sizeof(departments) / sizeof(departments[0]);
    int num_diagnoses = sizeof(diagnoses) / sizeof(diagnoses[0]);
    int num_treatments = sizeof(treatments) / sizeof(treatments[0]);

    // --- 文件操作和数据生成 ---
    FILE *fp = fopen("dataset.csv", "w");

    if (fp == NULL) {
        printf("Error: Unable to open the file.\n");
        return 1;
    }

    fprintf(fp, "record_id,t,k',doctor_name,patient_age,patient_gender,department,diagnosis,treatment\n");

    // 初始化随机数生成器
    srand(time(NULL));

    // 生成3000条随机数据
    for (int i = 1; i <= 3000; i++) {
        int record_id = i;
        int t = rand() % 7 + 1; 
        int k_prime = rand() % 6 + 1; 

        char doctor_name[50] = "";
        char patient_age[4] = "";
        char patient_gender[10] = "";
        char department[50] = "";
        char diagnosis[50] = "";
        char treatment[50] = "";

        // 创建一个索引数组，用于随机选择要填充的关键字
        int keyword_indices[] = {0, 1, 2, 3, 4, 5};
        
        // 打乱关键字索引数组以实现随机选择
        for (int j = 0; j < 6; j++) {
            int temp = keyword_indices[j];
            int random_index = rand() % 6;
            keyword_indices[j] = keyword_indices[random_index];
            keyword_indices[random_index] = temp;
        }
        
        // 根据 k' 的值，填充相应数量的关键字
        for (int j = 0; j < k_prime; j++) {
            switch (keyword_indices[j]) {
                case 0:
                    strcpy(doctor_name, doctor_names[rand() % num_doctors]);
                    break;
                case 1:
                    sprintf(patient_age, "%d", rand() % 80 + 1);
                    break;
                case 2:
                    strcpy(patient_gender, genders[rand() % num_genders]);
                    break;
                case 3:
                    strcpy(department, departments[rand() % num_departments]);
                    break;
                case 4:
                    strcpy(diagnosis, diagnoses[rand() % num_diagnoses]);
                    break;
                case 5:
                    strcpy(treatment, treatments[rand() % num_treatments]);
                    break;
            }
        }

        // 将生成的一行数据写入CSV文件
        fprintf(fp, "%d,%d,%d,%s,%s,%s,%s,%s,%s\n", 
                record_id, t, k_prime, doctor_name, patient_age, patient_gender,
                department, diagnosis, treatment);
    }

    fclose(fp);
    printf("dataset.csv file has been successfully generated with 1000 records from data pools.\n");

    return 0;
}
