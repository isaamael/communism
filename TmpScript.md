shell实现fa文件多行变一行

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' test.fa



sort排序时

-t $'\t'：指定TAB为分隔符
-k 1, 1: 按照第一列的值进行排序，如果只有一个1的话，相当于告诉sort从第一列开始直接到行尾排列
n:代表是数字顺序，默认情况下市字典序，如10<2