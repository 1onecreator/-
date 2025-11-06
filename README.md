本系统用c语言实现了在格密码背景下的针对于多关键词的有序级联搜索的实现
运行本系统的环境要求：Ubuntu20.04及以上版本
运行命令：
1.创建5000条数据数据
gcc -o CreateData CreateData.c
./Createdata
创建后的数据被保存在dataset.csv中
2.对5000条数据进行批量加密
gcc -o TestBigData TestBigData.c cipher.c generateKW.c global_sampling.c samplePre.c sampleLeft.c hash.c public_params.c -lm
./TestBigData
3.进行单条的数据加密或者搜索
gcc -o test2 test2.c Trapdoor.c verify.c cipher.c generateKW.c global_sampling.c samplePre.c sampleLeft.c hash.c public_params.c -lm -lpthread
./test2
