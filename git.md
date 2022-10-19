# Python基础语法学习

## 语法

- 注释：

  使用 #表示注释行,一般井号在行首，也可以在行的中间，那么#后面的为注释信息

  带有 #号的行python解释器会忽略

- 缩进

  使用缩进表示层次关系，或者区分不同的代码块，不像perl语言用{} 

  约定使用4个空格缩进

- 续行

  在行尾使用 \\

  如果使用各种括号，认为括号内是一个整体，括号内部跨行不用 \


Python示例代码：

```python
'''
Description: This script is used to onvert U in fasta sequence file to T
Date: 2019
Auther: omicsgene
'''

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import sys, os, argparse, os.path,re,math,time
#参数设置
parser = argparse.ArgumentParser(description='This script is \
    used to onvert U in fasta sequence file to T')
parser.add_argument('-f','--fasta',help='Please fasta file',required=True)
parser.add_argument('-o','--out_dir',help='Please input complete out_put directory path',
                    default = os.getcwd(),
                    required=False)

parser.add_argument('-n','--name',default ='demo_seq',required=False,
                    help='Please specify the output, demo_seq')
################################################################################
#读入参数，初始化路径
args = parser.parse_args()
dout=''
if os.path.exists(args.out_dir):
    dout=os.path.abspath(args.out_dir)
else:
    os.mkdir(args.out_dir)
    dout=os.path.abspath(args.out_dir)

output_handle = open(dout+'/'+args.name+'.fa', "w")
#循环处理序列将U转换成T,然后输出
for rec in SeqIO.parse(args.fasta, "fasta"):
    seq=rec.seq
    seq=str(seq.upper())
    seq=seq.replace("U", "T")
    seq_r = SeqRecord(Seq(seq,alphabet = IUPAC.IUPACAmbiguousDNA()), 
                      id=rec.id,
                      description=rec.description)
    SeqIO.write(seq_r, output_handle, "fasta")
output_handle.close()
```




## Python 标准数据类型和变量

#### 数字：

- 整数 integer：  11，223，-23   
- 布尔值 bool ：True  False
- 浮点数 float： 1.2  、 1.46 、  1.44e10 、 -1.6e-66
- 复数： 1+2j

#### 字符串 string：

- 使用 '  " 单引号或者双引号引用的字符
- ''' 和 """ 可以跨行，可以在其中任意的使用单双引号 

#### 空值

```
None
```


空值是Python里一个特殊的值，用None表示。None不能理解为0，因为0是有意义的，而None是一个特殊的空值。

#### 变量：


变量：

和初中代数的方程变量是一致的，只是在计算机程序中，变量不仅可以是数字，还可以是任意数据类型。

变量名：

1. 必须是大小写英文字母、数字和 _ 的组合，且不能用数字开头。例如，可将变量命名为message_1，但不能将其命名为1_message。
1. 变量名不能包含空格，但可使用下划线来分隔其中的单词。例如，变量名greeting_message可行，但变量名greeting message会引发错误。
1. 不要将Python关键字和函数名用作变量名，即不要使用Python保留用于特殊用途的单词，如print。
1. 变量名应既简短又具有描述性，让人一看就知道里面代表的是什么数据。例如，name比n好，student_name比s_n好，name_length比length_of_persons_name好。
1. 慎用小写字母l和大写字母O，因给他们可能被人错看成数字1和0；


```
a = 123 # a是整数
print(a)
a = 'ABC' # a变为字符串
print(a)
```

#### 变量的理解

```
a = 'abc'
b = a
a = 'xyz'
print(b)
```

1： 执行a = 'abc'，解释器创建了字符串'abc'和变量a，并把a指向'abc'：

```
graph LR
a-->abc
```

2：执行b = a，解释器创建了变量b，并把b指向a指向的字符串'abc'：

```
graph LR
a-->abc
b-->abc
```

3：执行a = 'xyz'，解释器创建了字符串'xyz'，并把a的指向改为'xyz'，但b并没有更改：

```
graph LR
a-->xyz
b-->abc
```


## Python 标准数据的操作


#### 数字操作运算

- 加减乘除  +  - *  /
- 求余运算 %    求幂运算  **
- 自然除/ 结果是浮点数，整除//。注意：在2.x中/和//都是整除。

注意：字符串可以相加和相乘

#### 字符串中特殊字符的转义：\

- 例如要在字符串中输入一个换行符： 

```
"hello python\n"
```

- 其他常见特殊意义的转义：

```
\\  \t  \r  \n  \'  \"
```

- 字符前缀 r  

```
print(r'my name is \n \"lucy\"')
```



练习：
打印以下变量

```
x = 122223
y = 45116e-789
a1 = 'Hello, world'
a2 = 'Hello, \'Lily\''
a3 = r'Hello, "Mike"'
a4 = r'''Hello,
omicsclass!'''
```

#### 字符串格式化输出

- %运算符格式化字符串
- %s表示用字符串替换，%d表示用整数替换，有几个%?占位符，后面就跟几个变量或者值，顺序要对应好。如果只有一个%?，括号可以省略。

```
print('Age: %s. Gender: %s' % (25, True))
print('growth rate: %d %%' % 7)
print('%2d-%02d' % (3, 1))
print('%.2f' % 3.1415926)
```



####  关于字符串操作的一些方法

==重点掌握：split，strip，join==

方法的调用用 . 操作符

| 方法                  | 说明                                           |
| --------------------- | ---------------------------------------------- |
| strip、rstrip、lstrip | 去除空白符。相当于对各个元素执行x.strip()      |
| split                 | 通过指定的分隔符将字符串拆分为一组子串         |
| join                  | 将字符串用作连接其他字符串序列的分隔符         |
| count                 | 返回子串在字符串中出现次数                     |
| endswith、startswith  | 如果字符串以某个后缀结尾（或开头），则返回True |
| lower、upper          | 分别将字母字符转换为小写或大写                 |
| replace               | 用另一个字符串替换指定子串                     |

   |  
ljust、rjust|	用空格(或其他字符)填充字符串的空白侧以返回符合最低宽度的字符串
index|	如果在字符串中找到子串，则返回子串第一个字符所在的位置。如果没有找到，则引发ValueError
find|	如果在字符串中找到子串，则返回第一个发现的子串第一个字符所在的位置。如果没有找到，则返回-1
rfind|	如果在字符串中找到子串，则返回最后一个发现的子串的第一个字符所在的位置。如果没有找到，则返回-1

示例代码：

```
s = 'A,T,C,G,A\n'
print(s)
print(s,end="")
print(s.strip())

s = 'A,T,C,G,A'
s.split(',')

"\t".join(["A","T","T"])


s = 'A,T,C,G,A'
s.replace('A','T')
s.count(',')
s.index('c')
s.find('c')
s.endwith("A")
s.lower()

#函数
len(s)


```

**注意:** 调用函数或者方法，一定要跟()


==**小知识**== 
Python中函数和方法的区别：

在Python中，**函数**(function) 和**方法**(method) 是有区别的：

- 所处的位置不同：函数是直接写在文件中而不是class中，方法是只能写在class中。

- 调用的方式不同：函数可以直接调用，例如：len()   print()等等，方法是在类中，一般用. 来调用，如： 前面学习的：lst.index() 



## python中的列表list和元组tuple

### 介绍

**list 列表**

- list是一种有序的集合，可以随时添加和删除其中的元素。
- 元素的类型可以是任意对象(数字，字符串，逻辑值，列表，对象等等)
- 使用[]表示列表


**tuple 元组**

- tuple：另一种有序列表,tuple和list非常类似，但是tuple一旦初始化就不能修改。
- 使用()表示元组

list的创建与初始化：

```
mylst=[]
mylst=[1,3,"a","xy",[11,22]]

```

tuple的创建：

```
classmates = ('Michael', 'Bob', 'Tracy')
t = ('a', 'b', ['A', 'B'])
```



### 列表list相关的操作&函数&方法

### 列表的查询：

##### 支持索引用中括号：[]

- 正索引：从左到右，从0开始
- 负索引：从右往左，从-1开始


```
lst=["A","T","A","a",22,[333,444]]

lst[1]
lst[0]
lst[-1]

```

##### 列表切片

```
lst=["A","T","A","a",22,[333,444]]
lst[0:3]
#表示，从索引0开始取，直到索引3为止，但不包括索引3。即索引0，1，2，正好是3个元素。
#如果第一个索引是0，还可以省略：
lst[:3]

lst[-2:-1]
lst[-2:]
```

##### 注意多维[][]

```
lst=["A","T","A","a",22,[333,444]]
lst[5][0]
```

### 列表元素修改

```
lst=["A","T","A","a",22,[333,444]]
lst[5][0]=333

lst[0]="omicsclass"

```


### 列表相关方法:

重点掌握方法：==append  extend  sort== 

| 方法                               | 功能                                                         |
| ---------------------------------- | ------------------------------------------------------------ |
| list.append(x)                     | 在列表末尾添加新的对象，只能添加一个元素                     |
| list.sort(key=None, reverse=False) | 对原列表进行排序                                             |
| list.extend(iterable)              | 在列表末尾一次性追加另一个序列中的多个值（用新列表扩展原来的列表） |
| list.count(x)                      | 统计某个元素在列表中出现的次数                               |
| list.index(x[, start[, end]])      | 从列表中找出某个值第一个匹配项的索引位置                     |
| list.insert(i, x)                  | 将对象插入列表                                               |
| list.pop([i])                      | 移除列表中的一个元素（默认最后一个元素），并且返回该元素的值 |
| list.remove(x)                     | 移除列表中某个值的第一个匹配项                               |
| list.reverse()                     | 反向列表中元素                                               |

### 示例代码

```
lst=["A","T","A","a",22,[333,444]]
#列表结尾追加元素
lst.append("omicsgene") 

#将可迭代的对象元素追加进来；
lst.extend(["omicsgene","omicsclass"]) 

#在指定的索引处插入元素
lst.insert(1,"omicsclass") 

#从左到右查询第一个匹配的元素删除
lst.remove("A")  

#不指定索引，就从列表结尾删除元素，并返回该元素
#也可以指定索引弹出；
lst.pop()
lst.pop(2)

#反转list，会修改原列表

lst.reverse()
```

### 列表相关的函数

| 函数           | 功能                   |
| -------------- | ---------------------- |
| len(list)      | 列表元素个数           |
| list(iterable) | 将可迭代对象转换为列表 |
| max(list)      | 返回列表元素最大值     |
| min(list)      | 返回列表元素最小值     |

```
a=[11,44,5,55,33,31,62,42]
len(a)
max(a)
min(b)

list([11,44,5,55,33,31,62,42])
list(range(10))

```

### 列表的排序

有全局函数sorted，和列表自带方法sort

**list.sort(key=None, reverse=False)**

- key 指定带有一个参数的函数或者方法，用于从每个列表元素中提取比较键(例如key=str.lower)。对应于列表中每一项的键会被计算一次，然后在整个排序过程中使用。默认值None表示直接对列表项排序而不计算一个单独的键值。
- reverse -- 排序规则，reverse = True 降序 ， reverse = False 升序（默认）。


```
## 1、最简单的排序
l = [5,2,3,1,4 ]
l.sort()
print(l)   

#反序
l.sort(reverse=True)
print(l)    

##2、字符串排序
StrList = ['Fast', 'Smooth', 'fast', 'isb', 'isa', 'smooth']
#一般字典序排列，但是大写在前，小写在后！！
StrList.sort()
print(StrList) 

##2.2忽略大小写，按abcd顺序
StrList.sort(key=str.lower)
print(StrList) 

##2.3按照字符串长度排序
StrList.sort(key=len)
print(StrList)

#一起使用两个参数
StrList.sort(key=len, reverse=True)
print(StrList) 
```


**sorted(iterable, key=None, reverse=False)**

- iterable -- 可迭代对象。 例如：列表  range对象  字符串等
- key -- 指定带有单个参数的函数，用于从 iterable 的每个元素中提取用于比较的键 (例如 key=str.lower)。 默认值为 None (直接比较元素)。
- reverse -- 排序规则，reverse = True 降序 ， reverse = False 升序（默认）。

```
a = [5,7,6,3,4,1,2]
b = sorted(a)       # 保留原列表
print(b)
b = sorted(a,reverse=True)  
print(b)
a = ['Fast', 'Smooth', 'fast', 'isb', 'isa', 'smooth']
b = sorted(a,key=str.lower,reverse=False)  
print(b)

```

**sort 与 sorted 区别：**

- 调用方式不同：
  sort 是应用在 list 上的方法，sorted 可以对所有可迭代的对象进行排序操作。

- 返回对象不同：
  list 的 sort 方法是在原来的列表上操作，无返回值（None），而函数 sorted 方法返回的是一个新的 list，而不是在原来的基础上进行的操作。



## python中的字典dict和集合set

- dict  字典

dict全称dictionary，使用键-值（key-value）存储，具有极快的查找速度。类似perl里面的hash。



### 字典的创建

字典的每个键值 key=>value对 用冒号 : 分割，每个键值对之间用逗号 , 分割，整个字典包括在花括号 {} 中 ,格式如下所示：


```
dict = {'a': 1, 'b': 2, 'b': '3'}
```

### 字典的特性：

1. 键必须不可变，所以可以用数字，字符串或元组充当，所以用列表就不行(列表是可变类型)：
1. 值可以没有限制地取任何python对象，既可以是标准的对象，也可以是用户定义的。


```
d = {'Name': 'Zara', 'Age': 7, 'Name': 'Manni'} 

#会报错
d = {['Name']: 'Zara', 'Age': 7}
```

### 修改字典

向字典中 添加新内容的方法是增加新的键/值对，修改或删除已有键/值对:


```
d = {'Name': 'Zara', 'Age': 7, 'Class': 'First'}
d['Name']                        #取值
d['Age'] = 18                    # 更新
d['School'] = "omicsclass"       # 添加
 
#取值打印
print( "d['Age']: ", d['Age'])
print("d['School']: ", d['School'])


# 删除字典元素


d = {'Name': 'Zara', 'Age': 7, 'Class': 'First'}
 
del d['Name']  # 删除键是'Name'的条目
d.clear()      # 清空字典所有条目
del d          # 删除字典
 
print("d['Age']: ", d['Age'])
print("d['School']: ", d['School'])
```

### 字典内置方法


| 方法                        | 描述                                          |
| --------------------------- | --------------------------------------------- |
| dict.keys()                 | 以列表返回一个字典所有的键                    |
| dict.values()               | 以列表返回字典中的所有值                      |
| dict.items()                | 以列表返回可遍历的(键, 值) 元组数组           |
| dict.get(key, default=None) | 返回指定键的值，如果值不在字典中返回default值 |

判断字典key是否存在

```

#使用格式：
#key in dict
#key not in dict

d = {'Name': 'Zara', 'Age': 7, 'Class': 'First'}
'Name' in d
'Name' not in d

```


### 集合

- set  集合

set 集合是一个无序的不重复元素序列。是一组key的集合，但不存储value。由于key不能重复，所以，在set中，没有重复的key。

- 集合创建


可以使用大括号 { } 或者 set() 函数创建集合。

注意：创建一个空集合必须用 set() 而不是 { }，因为 { } 是用来创建一个空字典。


#### 集合的创建与运算

```
genes = {'gene1', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5'}

genes = set(['gene1', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5'])

print(genes)   


#集合运算与判断
'gene1' in genes 

'gene11' in genes

#字符有列表特性，所以可以直接set   
a = set('omicsgene')
b = set('omicsclass')

# 集合a中包含而集合b中不包含的元素
a - b                              
# 集合a或b中包含的所有元素
a | b                              
# 集合a和b中都包含了的元素
a & b                              
# 不同时包含于a和b的元素
a ^ b                              

```


### 集合常用操作方法

| 方法           | 描述                                             |
| -------------- | ------------------------------------------------ |
| set.add(x)     | 参数作为一个元素添加到原set集合                  |
| set.update(s)  | 可以添加多个元素，且参数可以是列表，元组，字典等 |
| set.remove(x)  | 依据值删除，不存在会抛出异常KeyError             |
| set.discard(x) | 依据值删除，不存在不会报错                       |
| set.clear()    | 清空原set                                        |

代码示例：

```
basket = {'apple', 'orange', 'apple', 'pear', 'orange', 'banana'}
basket.add(("aa","dd"))
basket.update("dcc")
print(basket)


basket.remove("apple")
basket.remove("apple1")
basket.discard("apple1")

basket.clear()


```


## 数据类型及类型转换

Python3 支持 的标准数据类型总结

| 类型               | 类名                              |
| ------------------ | --------------------------------- |
| Number（数字）     | int、float、bool、complex（复数） |
| String（字符串）   | str                               |
| List（列表）       | list                              |
| Tuple（元组）      | tuple                             |
| Set（集合）        | set                               |
| Dictionary（字典） | dict                              |

Python3 的六个标准数据类型中：

- 不可变数据（3 个）：Number（数字）、String（字符串）、Tuple（元组）；
- 可变数据（3 个）：List（列表）、Dictionary（字典）、Set（集合）。

```
#批量赋值
a, b, c, d ,e= 20, 5.5, "ATGC", True, 4+3j

#获得数据类型
print(type(a), type(b), type(c), type(d), type(e))

#判断数据类型
a=1
isinstance(a, int)


```


类型转换函数：

| 函数           | 描述                    |
| -------------- | ----------------------- |
| int(x [,base]) | 将x转换为一个整数       |
| float(x)       | 将x转换到一个浮点数     |
| str(x)         | 将对象 x 转换为字符串   |
| list(s)        | 将序列 s 转换为一个列表 |
| set(s)         | 转换为可变集合          |
| tuple(s)       | 将序列 s 转换为一个元组 |

  |   
dict(d)	|创建一个字典。d 必须是一个 (key, value)元组序列。
complex(real [,imag])|	创建一个复数
repr(x)	|将对象 x 转换为表达式字符串
eval(str)|	用来计算在字符串中的有效Python表达式,并返回一个对象
frozenset(s)|	转换为不可变集合
chr(x)|	将一个整数转换为一个字符
ord(x)|	将一个字符转换为它的整数值
hex(x)	|将一个整数转换为一个十六进制字符串
oct(x)|	将一个整数转换为一个八进制字符串


## python中语句(if、for、while)


### if 语句 做判断

计算机之所以能做很多自动化的任务，因为它可以自己做条件判断。


#### python中常见比较判断运算符

==、  != 、 <  、> 、 >= 、 <=、in、not in

注意：in、not in ： 判断值是否存在列表中或者字典的键中


##### 示例代码：

```
#简单的判断
age = 20
if age >= 18:
    print('your age is', age)
    print('adult')

l=['Fast', 'Smooth', 'fast', 'isb', 'isa', 'smooth']
s = {'gene1', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5'}
d = {'gene1':122, 'gene2':2212, 'gene3':3121, 'gene4':2323, 'gene5':2543}


if 'Fast' not in l:
    print(True)

#else 语句
    
age = 3
if age >= 18:
    print('your age is', age)
    print('adult')
else:
    print('your age is', age)
    print('teenager')


#elif语句    
age = 3
if age >= 18:
    print('adult')
elif age >= 6:
    print('teenager')
else:
    print('kid')


#判断条件依次执行，有一个为真就跳出
age = 66
if age >= 18:
    print('adult')
elif age >= 6:
    print('teenager')
elif age >= 60:
    print('old people')
else:
    print('kid')
```

#### 逻辑运算符  增加判断条件

- 与 、或 、非 ：  and 、 or 、 not

1. and 如果第一个表达式为False，后面没有必要计算
1. or 如果第一个表达式是True，后面就没有必要计算了
1. 条件很多不知道判断优先顺序可以添加小括号


```
#逻辑运算符
age=11
gender="male"

if age>=18 and age<=60 and gender=="female":
    print("adult woman")

if age<18 or age >60 and gender=="female":
    print("kid or old woman")  
    
if age>=18 and age<=60 and not gender=="female":
    print("adult man")


```

#### 特殊数据真值表


| 对象/变量 | 值   |
| --------- | ---- |
| ""        | 假   |
| "string"  | 真   |
| 0         | 假   |
| () 空元组 | 假   |
| [] 空列表 | 假   |
| {} 空字典 | 假   |
| None      | 假   |

### for 循环语句

计算机之所以能批量的执行数据分析，就是因为有循环语句。


for...in循环，in 后面跟可迭代对象


```
#列表遍历
names = ['Michael', 'Bob', 'Tracy']
for name in names:
    print(name)
    
    
sum = 0
for x in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    sum = sum + x
print(sum)    


#产生连续的数字函数：range(start, stop[, step])

range(10)
range(2,10)
#也可以转换成list对象
list(range(10))
list(range(1, 11))
list(range(0, 30, 5))
list(range(0, 10, 3))
list(range(0, -10, -1))
list(range(0))
list(range(1, 0))


sum = 0
for x in range(101):
    sum = sum + x
print(sum)


sum = 0
for x in range(3,101):
    sum = sum + x
print(sum)


sum = 0
for x in range(3,101,2):
    sum = sum + x
print(sum)


#enumerate函数  同时获取索引还有值

#用法enumerate(sequence, [start=0])

names = ['Michael', 'Bob', 'Tracy']
for index,value in enumerate(names):
    print(index,value)


#字典遍历
dict = {'Name': 'Zara', 'Age': 7, 'Class': 'First'}
for k,v in dict.items():
    print(k,v)

dict = {'Name': 'Zara', 'Age': 7, 'Class': 'First'}
for k in dict.keys():
    print(k,dict[k])
```

### while 循环语句

while循环，只要条件满足，就不断循环，条件不满足时退出循环


```
sum = 0
n = 99
while n > 0:
    sum = sum + n
    n = n - 2
print(sum)
```

### 循环控制语句

break  提前退出循环
continue  提前结束本轮循环



```
n = 1
while n <= 100:
    if n > 10: # 当n = 11时，条件满足，执行break语句
        break # break语句会结束当前循环
    print(n)
    n = n + 1
print('END')


n = 0
while n < 10:
    n = n + 1
    if n % 2 == 0: # 如果n是偶数，执行continue语句
        continue # continue语句会直接继续下一轮循环，后续的print()语句不会执行
    print(n)
```


==练习题==

1. for 循环打印星号金字塔：



```
     *
    ***
   *****
  *******
 *********
```

答案：



```
for i in range(1, 6):
     for j in range(0, 6 - i):
         print (" ",end="")
     
     print( "*"*(i),end="")
     print( "*"*(i-1),end="")       
     print("")
```

2. 利用for循环，打印乘法口诀：

```
1*1=1	
1*2=2	2*2=4	
1*3=3	2*3=6	3*3=9	
1*4=4	2*4=8	3*4=12	4*4=16	
1*5=5	2*5=10	3*5=15	4*5=20	5*5=25	
1*6=6	2*6=12	3*6=18	4*6=24	5*6=30	6*6=36	
1*7=7	2*7=14	3*7=21	4*7=28	5*7=35	6*7=42	7*7=49	
1*8=8	2*8=16	3*8=24	4*8=32	5*8=40	6*8=48	7*8=56	8*8=64	
1*9=9	2*9=18	3*9=27	4*9=36	5*9=45	6*9=54	7*9=63	8*9=72	9*9=81	
```

答案：

```
for m in range(1, 10):
    for n in range(1, m+1):
        print("%d*%d=%d\t"%(n,m,n*m), end="")
    print("")
```


## python 读写数据文件

### python文件读写分为三步

1. 打开文件获取文件对象
1. 操作文件
1. 关闭文件

```
fr = open("test.txt","r")      #打开文件
ff = fr.read()                 #读取文件所有内容  （不建议使用，如果文件内容巨大，内存会爆）
print(ff) 
```

### 文件类型

| 参数 | 描述                                                         |
| ---- | ------------------------------------------------------------ |
| r    | 只读，默认模式   打开数据文件                                |
| w    | 只写，不可读，若文件不存在则创建，若存在，则删除内容，写入新内容 |
| a    | 只追加，不可读，若文件不存在则创建，存在则追加新内容         |


### 文件方法（不常用）

 前三个方法：在大文件时慎用，会把内容读到内存中，占用大内存

| 方法      | 描述                                                         |
| --------- | ------------------------------------------------------------ |
| f.write() | 字符串写入一个打开的文件                                     |
| f.close() | 刷新缓冲区里任何还没写入的信息，并关闭该文件，这之后便不能再进行写入。 |

|
f.read()       |    #读取所有内容
f.readline()    |   #读取一行
f.readlines()  |    #读取所有文件内容，返回一个list
f.seek(0)    |      #当前文件指针位置在0位
f.writelines(["a","b"])  |  #把列表写入文件




#### 文件的读入与写出（常用方法）：

文件的读入与写出主要用的是for循环：

```
fr=open("input.fa","r")   #输入文件，读取数据
fw=open("output.fa","w")   #输出文件，写出结果
for line in fr:              #循环一行一行的读取文件                                 
    new_line = line.replace("A","T")                         
    fw.write(new_line)
fr.close()  #关闭文件
fw.close()  #关闭文件
```


测试文件内容：
input.fa

```
AGTTAGCGGATAATGGCCATCAAAGCAACGCTTACCAACACTGCACCCCTTGTTTTGGAAATGCAACCAC
AAAAGCATTGGACACTTGCTTACTTCAAATAAAACACATTTAAACAATTAGATGACGTGATGGACCAGAA
TGGCGCATCGGGAAGTCATCCGAACAGGCTATCCCAAGGAAGAGGAGCCCATGCGCGCGAACGTGGCGCC
ACAGTTTCCGCGGCGGCAAATCGGAGTAACATTATCGACGAAATGGCCAAAATATGCGAAGCCGATCGCC
AGACTTTCGCCATCGCTCGACGGACTCGGGGTCACGAGCGGCTTGCGGTGGACAACAGCGACTTCGTCGC
CGTGGAGGATCTTATTTTGTCCTACGCAGAGCCCACGCCCGAGGACCAGGTCGAGATGATCATGAGCGAC
TTTTGCTCGTCTCCAACATACGCAGAGGATGAGGATGAGCCCAGCCATGAGTCGGAGCCGTGGTTTCGAT
TTCGCAACAAAAGGATCAGAACCTACAGCCGGAAGAGGGATCCCAAAAGCCACAAGGCCGTTCAAAACGA
GAAGCGTAGAGGTTCCTCAGGCCTCTCCGTGCAGAGGGATCTCAATACTTCGTTCACATCTATGGCTTGT
GATTTCGATGCTTCATCACAGAAGATACACGAGGTCCTTTTGAACCTCAGTCAATACTTTTCCGCGACCG
CGACAGCTTCCGGTCCGACTCCTGTCCCATCGCAAATAGATCTGCCAACCGAAGCAAGGCAGGATT
```

---

==学习任务==

- 编写脚本要求：
- [x] 读取gff文件，筛选出文件中1号染色体100000-500000之间的基因
- [x] 输出基因的名字，染色体，起始位置，终止位置信息用tab分隔各列

输入文件下载地址：

ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.chromosome.1.gff3.gz


注意：几乎所有的生物数据文件都是文本文件，都可以用notepad++，editplus等文本编辑器打开，不要被文件后缀名迷惑，我们的编程语言python绝大多数情况下打开处理的也都是文本文件。



---

**答案：**

```
fr=open("D:\python_script\Arabidopsis_thaliana.TAIR10.44.chromosome.1.gff3\Arabidopsis_thaliana.TAIR10.44.chromosome.1.gff3","r")
fw=open("D:\python_script\Arabidopsis_thaliana.TAIR10.44.chromosome.1.gff3\out.txt","w")

for line in fr:
    line=line.strip()
    if not line[0]=="#":
        tmp=line.split("\t")
        if tmp[0]=="1" and tmp[2]=="gene" and int(tmp[3])>100000 and int(tmp[3])<500000:
            geneID=tmp[8].split(";")[0].split("=")[1]
            #mystr=tmp[0]+"\t"+tmp[3]+"\t"+tmp[4]+"\t"+geneID+"\n"
            mystr="\t".join([tmp[0],tmp[3],tmp[4],geneID])+"\n"
            fw.write(mystr)
            #fw.write(line+"\n")
    
fr.close()
fw.close()
```