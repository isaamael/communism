

# Python高级编程

## 1 函数

### 1.1 Python内置函数：

Python内置了很多有用的函数，我们可以直接调用。

内置函数总结：https://docs.Python.org/zh-cn/3.7/library/functions.html

```
abs(-100)
max(1, 2)
int('123')
float('12.34')
sum([2,323,23])
```

### 1.2 Python 用户自定义函数

自定义函数：把具有独立功能的代码块组织成为一个小模块。

好处：
1. 代码复用，提高编程效率，使程序易于理解。
2. 自己学会定义函数，有助于理解看懂别人函数的帮助
#### 定义一个函数
你可以定义一个自定义功能的函数，以下是简单的规则：

1.  函数代码块以 def 关键词开头，后接函数标识符名称和圆括号()。
1.  括号里面可设置参数，参数不是必必需的(可有可无)，传入参数必须放在圆括号中间。
1.  函数的第一行语句可以选择性地（可有可无）使用文档字符串—用于存放函数使用说明。
1.  函数内容以冒号起始，并且缩进。
1.  return [表达式] 结束函数，选择性地返回一个值给调用方。不带表达式的return相当于返回 None。

参数作用：增加函数的通用性，针对相同的数据处理逻辑，能够适应更多的数据
1. 在函数内部，把参数当作变量使用，进行需要的数据处理
2. 函数调用时，按照函数定义的参数顺序，把希望在函数内部处理的数据，通过参数传递

#### 示例

```

#定义
def hello(s):
    '''this function is used to say hello,
    and exit.
    
    '''
    print('hello '+s)
    return
    
#调用
hello(s="omicsgene")

```


#### 不同类型参数的设置与使用

##### 参数根据位置传递

如果我们设置了多个参数，python解释器会根据位置进行传递；

```

def power(x, n):
    s = 1
    while n > 0:
        n = n - 1
        s = s * x
    return s

power(2,3)

power(x=2,n=3)

```


##### 必需参数

直接设定参数，没有给参数的默认值，函数调用时必需给参数传入值

```

def power(x, n):
    s = 1
    while n > 0:
        n = n - 1
        s = s * x
    return s
    
power(4)
```



##### 设置默认参数(函数调用时可有可无)
设置参数是可以设定参数的默认值，函数调用时可选给参数传入值们不给就是默认值。

```
def power(x, n=2):
    s = 1
    while n > 0:
        n = n - 1
        s = s * x
    return s

def info(name, gender, age=6, city='Beijing'):
    print('name:', name)
    print('gender:', gender)
    print('age:', age)
    print('city:', city)
```


##### 可变参数位置参数 

\* 在参数名字前面加一个星号，使参数成为一个可变参数，参数传入一个列表

```
def add(a,b):
    return a+b

def add1(*numbers):
    sum = 0
    for n in numbers:
        sum = sum + n 
    return sum
```


##### 可变关键字参数 

** 在参数名字前面加两个星号，使参数成为一个可变关键字参数，参数传入一个字典
```
def person(name, age, **kwargs):
    print('name:', name, 'age:', age, 'other:', kwargs)


#联合使用，参数如何传递？  
def foo(a,b,c=22,*args,**kwargs):
    print("a:",a,",","b:",b,",","c:",c,",","args:",args,",","kwargs:",kwargs)
```

#####  强制  key=value形式传递

有个参数只有一个星号,表示星号后面的参数必须用（key=value）的形式进行传递参数，不能省略key。

```
def foo1(a,b,*,c=3):
    print(a,b,c)
    
def foo2(a,b,c=3):
    print(a,b,c)

foo1(1,2,4)
foo2(1,2,4)
```


##### 函数编写与调用注意事项

1. 必需参数一般写在前面，有默认值的参数紧随其后（位置参数）
1. 参数名省略时，注意位置

```
def foo(a,b,c=22,*args,**kwargs):
    print(a,b,c,args,kwargs)

foo()
```



### 1.3 Python中匿名函数


有时候一个函数只是临时使用一下，以后就不会再使用了，要是用def定义比较麻烦。这个时候就可以用lambda来定义一个匿名函数

##### 语法：

变量名= lambda [arg1[, arg2, ... argN]]: expression
##### 注意:

1. 参数：可选，通常以逗号分隔的变量表达式形式，也就是位置参数
1. 表达式中不能包含 循环，return
2. 可以包含 if...else...语句.
3. 表达式计算的结果直接返回


```

# 使用lambda的表达式
lambda x, y: x + y

# 使用def定义的函数
def add( x, y ):
  return x + y

# lambda也允许有默认值和使用变长参数
lambda x, y = 2: x + y
lambda *z: z

# 调用lambda函数
a = lambda x, y: x + y
a(1, 3)

b = lambda x, y = 2: x + y
b(1)

b(1, 3)

c = lambda *z: z
c(10, 'test')


```



#### 匿名函数应用


```
#自定义排序关键字

students = [('john', 'A', 15), ('jane', 'B', 12), ('dave', 'B', 10)]
sorted(students, key=lambda s: s[2])            # 按年龄排序
sorted(students, key=lambda s: s[2], reverse=True)       # 按降序


#字典排序：


dic = {'a':2,'b':1}
#按照key排序
d = sorted(dic.items(), key = lambda k:k[0])

#按照values排序
 e = sorted(dic.items(), key = lambda k:k[1])


```
### 1.4 变量作用域
变量的作用域决定了在哪一部分程序你可以访问那个特定的变量名称。一个程序的所有的变量并不是在任何位置都可以访问的。访问权限决定于这个变量是在哪里赋值的。

两种最基本的变量作用域如下：

- 全局变量
- 局部变量


##### 全局变量和局部变量

- 定义在函数内部的变量拥有一个局部作用域，定义在函数外的拥有全局作用域。

- 局部变量只能在其被声明的函数内部访问，而全局变量可以在整个程序范围内访问。



```
total = 0; # 这是一个全局变量

def mysum( arg1, arg2 ):
   total = arg1 + arg2; # total在这里是局部变量.
   print("inside  ", total)
   return total
 
#调用mysum函数
mysum( 10, 20 )
print("global var ", total)
```


​    
### 1.5 变量传递给参数参数的可变性

函数的参数传递中传入可变数据和不可变数据会有不同：

- 不可变数据（3 个）：Number（数字）、String（字符串）、Tuple（元组）；
- 可变数据（3 个）：List（列表）、Dictionary（字典）、Set（集合）。

##### 传不可变对象实例,传进去的变量不改变原变量

```
def ChangeInt(a):
    a = 10
 
b = 2
ChangeInt(b)
print (b) # 结果是 2
```

##### 传可变对象实例,传进去的变量可改变原变量

```
def changeme(mylist):
   "修改传入的列表"
   mylist.append([1,2,3,4]);
   print("Inside value:", mylist)
   return
 
# 调用changeme函数
mynum = [10,20,30]
changeme(mynum)
print("Outside value:", mynum)
```

---

## 2 面向对象的编程与获取帮助
### 2.1 类和实例

>面向对象最重要的概念就是类（Class）和实例（Instance），类描述了一组具有相同特性（属性）和相同行为（方法）的对象(Object)。面向对象的编程语言最大的特色就是可以编写自己所需的数据类型，以更好的解决问题。

>必须牢记类是抽象的模板，而实例是根据类创建出来的一个个具体的“对象”，每个对象都拥有相同的方法或者属性，但各自的数据可能不同。

![image](https://note.youdao.com/yws/api/group/92990402/noteresource/8DB0B3FB6790404EADCA0B71D9EB70F1/version/987?method=get-resource&shareToken=08DB9C16FB6D47AA8EBD104017423D2D&entryId=446741437)

#### 面向对象编程特点

面向对象编程(Object-oriented programming，缩写：OOP)的3个基本特征是：封装、继承、多态

- 封装：将属性和方法(数据和功能)封装在一起形成类。 

- 继承：可以使用现有类的功能，并在无需重新编写原来的类的情况下对这些功能进行扩展。

- 多态：允许让父类的指针分别指向不同的子类, 调用不同子类的同一个方法, 会有不同的执行效果



### 2.2 面向对象编程

#### 定义一个类

##### 类的组成

- 属性(对象的属性) ——变量：状态、静态的
- 方法(对象的行为) ——函数：过程、动态的

##### 类的方法与属性


>在类的内部，使用 def 关键字来定义一个方法，类的方法与普通的函数只有一个特别的区别——他们的第一个参数必须是 self。

>self代表类的实例，而非类。

##### 构造方法

类有一个名为 \_\_init__() 的特殊方法（构造方法），该方法在类实例化时会自动调用。


```
#类的定义
class Car:
    '''this class define a car '''
    #类属性  共有属性
    wheels=4
    #构造方法
    def __init__(self, make, model, year):
        #成员属性
        self.make = make
        self.model = model
        self.year = year
        self.orometer_reading = 0
    #类方法
    def get_description(self):
        long_name = str(self.year) + ' ' + self.make + ' ' + self.model+" "+str(self.wheels)
        return long_name
 
    def get_odometer(self):
        print("This car has "+ str(self.orometer_reading) + " miles on it")
 
    def increase(self,miles):
        self.orometer_reading +=miles

```

#### 创建类的实例化对象，访问属性和使用方法
现在让我们新建一个对象my_car:
```
#类实例my_car
my_car = Car("yellow", "beetle", 1967)
#查看属性
print(f" My {my_car.color} car {my_car.model} is made in {my_car.year}")

#属性修改
my_car.color="black"

#对象方法调用

my_car.get_description()
my_car.increase(1000)
my_car.read_odometer()
my_car.update_orometer(20000)

#查看类或者实例所有的属性与方法
dir(my_car)
dir(Car)
```
### 2.3 Python 类定义中的特殊属性与方法

Python中用下划线作为变量前缀和后缀指定特殊变量


- _xxx 不能用’from module import *’导入
- \_\_xxx 类中的私有变量名
- \_\_xxx__ 系统定义的名字
- 核心风格：避免用下划线作为变量名的开始。


因为下划线对解释器有特殊的意义，而且是内建标识符所使用的符号，我们建议程序员避免用下划线作为变量名的开始。

==一般来讲，变量名_xxx被看作是“私有的”，在模块或类外不可以使用。==

当变量是私有的时候，用_xxx 来表示变量是很好的习惯。因为变量名__xxx__对Python 来说有特殊含义，对于普通的变量应当避免这种命名风格。

- “单下划线” 开始的成员变量叫做保护变量，意思是只有类对象和子类对象自己能访问到这些变量；
- “双下划线” 开始的是私有成员，意思是只有类对象自己能访问，连子类对象也不能访问到这个数据。

##### Python类的特殊属性
属性|	含义
--|--
\_\_class__|对象或类所属的类
\_\_name__|类、函数、方法等的名字
\_\_dict__|类的属性  以key-value形式展示的字典，展示出属性所对应的值
\_\_module__|类定义所在的模块名称
\_\_doc__|类、函数的文档字符串，如果没有定义则为None

```
#print(my_car.__name__)
print(my_car.__doc__)
print(my_car.__dict__)
print(my_car.__module__)
print(my_car.__class__)


print(Car.__name__)
print(Car.__doc__)
print(Car.__dict__)
print(Car.__module__)
print(Car.__class__)
```



##### Python类特殊方法
方法|	功能说明
--|--
\_\_new__()	|类的静态方法，用于确定是否创建对象
\_\_init__()|	构造函数，生成对象时调用
\_\_dir__|返回类或者对象的所有方法与属性，dir()操作实例就是调用
\_\_del__()|	析构函数，释放对象时调用
\_\_add__()|	+
\_\_sub__()|	-
\_\_mul__()|	*
\_\_truediv__()|	/
\_\_floordiv__()|	//
\_\_mod__()|	%
\_\_pow__()|	**
\_\_repr__()	|打印,转换
\_\_setitem__()|	按照索引赋值
\_\_getitem__()|	按照索引获取值
\_\_len__()|	计算长度
\_\_call__()|	函数调用
\_\_contains__() |	in
\_\_eq__()|	==
\_\_ne__()|	!=
\_\_lt__()|	<
\_\_le__()|	<==
\_\_gt__()|	>
\_\_ge__()|>=
\_\_str__()|	转换为字符串
\_\_shift__(), \_\_rshift__()|	<<, >>
\_\_and__(), \_\_or__()|	&, |
\_\_invert__(), \_\_xor__()|	~, ^
\_\_iadd__(), \_\_isub__()|	+=, -=


##### 例子：实现特殊方法__repr__：


```

## __repr__
class Person:
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return "hello %s." % self.name

 p = Person('hkey')

 p
 print(p)

```

定义了__repr__方法，不管是直接打印对象还是通过print打印对象，都是走的__repr__中定义的格式。


更多类的特殊方法了解：https://www.omicsclass.com/article/1033

### 2.4 类的继承 (选修)

```
#类定义
class people:
    #定义属性
    name = ''
    age = 0
    #定义私有属性,私有属性在类外部无法直接进行访问
    __weight = 0
    #定义构造方法
    def __init__(self,n,a,w):
        self.name = n
        self.age = a
        self.__weight = w
    #定义类方法
    def speak(self):
        print("%s speak: I am %d years old." %(self.name,self.age))
 
#单继承示例
class student(people):
    grade = ''
    def __init__(self,n,a,w,g):
        #调用父类的构函
        people.__init__(self,n,a,w)
        self.grade = g
    #覆写父类的方法
    def speak(self):
        print("%s speak: I am %d years old， I am in %d grade of primary school."%(self.name,self.age,self.grade))
 
 
 
s = student('ken',10,60,3)
s.speak()
```
### 2.5 Python获取帮助

#### Python 官方中文帮助
https://docs.Python.org/zh-cn/3.7/

#### 编辑器提供帮助
pycharm中的documentation内置显示（默认快捷键为Ctrl+Q），选中函数，Ctrl+Q如下：

也有External documetation，快捷键为Shift+F1

#### 陌生类获取类的属性和方法

dir()函数主要用来查看对象的属性和方法，再去了解这些属性和方法，学习类的使用。

####  __doc__属性查看帮助文档


任务答案




```
a = [1,2,3]
a.reverse.__doc__
```


==任务答案：==

```
def my_abs(x):
    a=0
    if x>=0:
        a=x
    else:
        a=0-x
    return a
    
my_abs(-100)
```


---

## 3 Python 标准库中包的使用

python标准库帮助地址：https://docs.python.org/zh-cn/3.7/library/index.html


### 3.1 包常见导入方法

```

# 1 直接导入
import os
import time
import sys
import os,time,sys,re  #每个包之间用逗号隔开；

#导入的同时改名字
import sys as system


#2 导入指定的方法，模块等

#导入包里面的指定的函数或者类
from os import path
from os import path, walk, unlinkfrom

#导入包里面所有的内容
from os import *

```

这些导入的包，其实就是别人写好的代码文件(*.py)，我们导入到我们的程序中就可以直接使用里面的方法，函数，类等，省去自己编写代码的麻烦。如果了解Python包导入的机制，自己也可以写一些包，共享给别人，提供代码的复用性，从而提高我们的开发效率。

### 3.2 Python 自定义模块导入讲解

#### 模块搜索路径
导入过程首先需要定位导入文件的位置，也就是，告诉Python到何处去找到要导入的文件，因此，需要设置模块的搜索路径。在大多数情况下，Python会自动到默认的目录下去搜索模块；如果要在默认的目录之外导入模块，就需要知道Pyhon搜索模块路径的机制。

Python搜索模块的路径是由四部分构成的：==程序的主目录、PATHONPATH目录、标准链接库目录和.pth文件的目录，这四部分的路径都存储在sys.path 列表中。==

##### 1， 程序的主目录

主目录是指程序所在的目录，Python首先会到主目录中搜索模块。

因为主目录总是第一个被搜索，如果模块完全处于主目录中，所有的导入都会自动完成，而不需要单独配置路径。

##### 2，PATHONPATH目录

PythonPATH目录是指PythonPATH环境变量中配置的目录，是第二个被搜索的目录，Python会从左到右搜索PythonPATH环境变量中设置的所有目录。

##### 3，标准链接库目录

标准链接库目录是Python按照标准模块的目录，是在安装Python时自动创建的目录，通常不需要添加到PythonPATH目录中。例如：

##### 4，路径文件（.pth文件）

在模块搜索目录中，创建路径文件，后缀名为.pth，该文件每一行都是一个有效的目录。Python会读取路径文件中的内容，每行都作为一个有效的目录，加载到模块搜索路径列表中。简而言之，当路径文件存放到搜索路径中时，其作用和PYT)HONPATH环境变量的作用相同。


### 3.3 导入自定义模块指定路径的方法

##### 方法1：简单不用操作
把自己写的模块，复制到自己的程序同级目录就可以直接导入；


##### 方法2: 函数添加  （临时添加）
1. import sys
2. 查看sys.path
3. 添加sys.path.append("D:\\\perl_script")

##### 方法3: 增加.pth文件，推荐！   （永久添加）
- 在site-packages目录添加一个文件如mypkpath.pth，必须以.pth为后缀，文件内写上你要加入的模块文件所在的路径，可以添加多行。

- 如果运行在Windows和Python3.0中，如果Python安装目录的顶层是C:\Users\Administrator\AppData\Local\Programs\Python\Python37，那么可以把自定义的路径文件 mypath.pth 放到该目录中。

- 也可以放到标准库所在位置的site-packages子目录中（C:\Users\Administrator\AppData\Local\Programs\Python\Python37\Lib\site-packages），来扩展模块的搜搜路径。

##### 方法4: 修改环境变量  （永久添加）
- windows用户可以修改系统环境变量 PYTHONPATH
- linux用户也可以添加环境变量 PYTHONPATH


### 练习：

导入自己写的一个模块：

把以下内容存储到文件cars.py中，然后用import导入该模块：

```

CAR_num=11

make=["Ford","Rolls-royce","Volkswagen"]

def mysum( arg1, arg2 ):
   total = arg1 + arg2; # total在这里是局部变量.
   print("inside  ", total)
   return total;


def das_auto(make="Volkswagen", model="Magotan", year=2019):
    '''this function is used to create Volkswagen car'''
    my_car=Car("Volkswagen","Magotan",2019)
    return my_car


#类的定义
class Car:
    '''this class define a car '''
    #类属性  共有属性
    wheels=4
    #构造方法
    def __init__(self, make, model, year):
        #成员属性
        self.make = make
        self.model = model
        self.year = year
        self.orometer_reading = 0
    #类方法
    def get_description(self):
        long_name = str(self.year) + ' ' + self.make + ' ' + self.model+" "+str(self.wheels)
        return long_name
 
    def read_odometer(self):
        print("This car has "+ str(self.orometer_reading) + " miles on it")
 
    def update_orometer(self,miles):
        if miles >= self.orometer_reading:
            self.orometer_reading = miles
        else:
            print("You can'troll back an odometer")
 
    def increase(self,miles):
        self.orometer_reading +=miles

```

##### ==小知识== 模块和包的区别

（1）模块：是一个单独的.py文件，用于存放一些功能相关的代码，可以是代码更加容易维护，提高代码的重用价值

（2）包：是一个有层级的目录结构，包含n个模块或者n个子包，包中一定要有__init__.py文件


### 3.4 Python中 os sys time math包的学习


os 常用方法


方法|	 说明
--|--
os.getcwd()|得到当前工作目录，即当前Python脚本工作的目录路径。
os.mkdir() |方法用于以数字权限模式创建目录，
os.listdir(path)|返回指定目录下的所有文件和目录名。
os.walk()| 方法用于通过在目录树中游走输出在目录中的文件名，向上或者向下。
os.remove(path)|方法用来删除一个文件。如果指定的路径是一个目录，将抛出OSError，在Unix, Windows中有效
os.rmdir(path)|方法用于删除指定路径的目录。仅当这文件夹是空的才可以, 否则, 抛出OSError。
os.system(command)|函数用来运行shell命令。
os.linesep|字符串给出当前平台使用的行终止符。例如，Windows使用'\r\n'，Linux使用'\n'而Mac使用'\r'。
os.sep|可以取代操作系统特定的路径分隔符。windows下为 “\\\”,linux 下为"/"
os.chdir(dirname)|改变工作目录到dirname
os.path.isfile(path)|方法分别检验给出的路径是否一个文件，输入路径必须是绝对路径。
os.path.isdir(path)|方法检验给出的路径是否一个目录，输入路径必须绝对路径。
os.path.exists()|方法用来检验给出的路径是否真地存在
os.path.join(path,name)|连接目录与文件名或目录;使用“\”连接
os.path.basename(path)|返回文件名
os.path.dirname(path)|返回文件路径
os.path.abspath(name)|获得绝对路径
os.path.getsize(name)|获得文件大小
os.path.split(path) |将path分割成目录和文件名二元组返回。
os.path.splitext()|分离文件名与扩展名


#### 包使用示例代码

```
#遍历一个目录，获得所有文件的路径

import os  
os.chdir("D:\Python_script")
cwd = os.getcwd()  
for dir_path, dir_names, file_names in os.walk(cwd): 
    for file_name in file_names:    
        p=os.path.join(dir_path,file_name)
        print(p)
    for dir_name in dir_names:
        p=os.path.join(dir_path,dir_name)
        print(p)


#windows系统中调用命令
os.system('copy a.v b.v')
os.system('del b.v')
os.system('rename a.v b.v')


#linux系统中调用命令
os.system('cp a.v b.v')
os.system('rm b.v')
os.system('mv a.v b.v')

```


#### sys包常用属性或者方法
方法|	 说明
--|--
sys.argv |命令行参数List，第一个元素是程序本身路径
sys.exit(n) |退出程序，正常退出时exit(0)
sys.version |获取Python解释程序的版本信息
sys.modules |返回系统导入的模块字段，key是模块名，value是模块
sys.path |返回模块的搜索路径，初始化时使用PythonPATH环境变量的值
sys.stdout |标准输出
sys.stdin |标准输入
sys.stderr |错误输出

#### math 常用属性和方法
方法|	 说明|示例
--|--|--
math.e	 |自然常数e	| >>> math.e
math.pi	 |圆周率pi	 |>>> math.pi
math.log10(x)	 |返回x的以10为底的对数	 |>>> math.log10(2)
math.pow(x, y)	 |返回x的y次方|>>> math.pow(5,3)
math.sqrt(x)	 |返回x的平方根	 |>>> math.sqrt(3)
math.ceil(x)	 |返回不小于x的整数	 |>>> math.ceil(5.2)
math.floor(x)	 |返回不大于x的整数	 |>>> math.floor(5.8)
math.fabs(x)	| 返回x的绝对值	 |>>> math.fabs(-5)




---
练习：
将文件打开然后，批量的合并到一个文件中

```
fw=open("all.fa","w")
for i in os.listdir(os.getcwd()):
    if i.endswith("fa"):
        f=open(i,"r")
        for line in f:
            fw.write(line)
        f.close()
fw.close() 
```


## 4 正则表达式：
### 4.1 了解正则表达式

正则表达式是用于处理字符串（String）的强大工具，拥有自己独特的语法以及一个独立的处理引擎，效率上可能不如str自带的方法，但功能十分强大。得益于这一点，在提供了正则表达式的语言里，正则表达式的语法都是一样的，区别只在于不同的编程语言实现支持的语法数量不同；

#### 正则表达式功能
- 查找/匹配               
- 替换                          
- 捕获
- 计数



#### 正则表达式通配符
- 通配符是特殊字符，在正则表达式中有特殊意义
- 直接的单词或者数字是原样匹配


#### 常见通配符学习
模式匹配通配符|	描述
--|--
^	|匹配字符串的开头
$|	匹配字符串的末尾。
.|	匹配任意字符，除了换行符，当re.DOTALL标记被指定时，则可以匹配包括换行符的任意字符。
[...]|	用来表示一组字符,单独列出：[amk] 匹配 'a'，'m'或'k'，[a-zA-Z0-9]
[^...]	|不在[]中的字符：[^abc] 匹配除了a,b,c之外的字符。
\w	|匹配字母数字及下划线
\W	|匹配非字母数字及下划线
\s	|匹配任意空白字符，等价于 [\t\n\r\f].
\S	|匹配任意非空字符
\d	|匹配任意数字，等价于 [0-9].
\D	|匹配任意非数字
*	|匹配0个或多个的表达式,如:be*:  b, be,beeeee, bee。
+	|匹配1个或多个的表达式,如: be+ : be,beeeee, bee。
?	|匹配0个或1个由前面的正则表达式定义的片段，非贪婪方式 .如  bo? 只能匹配：b,bo
{ n}	|精确匹配 n 个前面表达式。例如， o{2} 不能匹配 "Bob" 中的 "o"，但是能匹配 "food" 中的两个 o。
{ n,}	|匹配 n 个前面表达式。例如， o{2,} 不能匹配"Bob"中的"o"，但能匹配 "foooood"中的所有 o。"o{1,}" 等价于 "o+"。"o{0,}" 则等价于 "o*"。
{ n, m}|	匹配 n 到 m 次由前面的正则表达式定义的片段，贪婪方式
\| 	|或者的意思，如：a\| b  匹配a或b
()	|对正则表达式分组并记住匹配的文本,正则表达式捕获功能




##### 选修内容
模式匹配通配符|	描述
--|--
(?imx)	|正则表达式包含三种可选标志：i, m, 或 x 。只影响括号中的区域。
(?-imx)|	正则表达式关闭 i, m, 或 x 可选标志。只影响括号中的区域。
(?: re)	|类似 (...), 但是不表示一个组
(?imx: re)	|在括号中使用i, m, 或 x 可选标志
(?-imx: re)	|在括号中不使用i, m, 或 x 可选标志
(?#...)	|注释.
(?= re)	|前向肯定界定符。如果所含正则表达式，以 ... 表示，在当前位置成功匹配时成功，否则失败。但一旦所含表达式已经尝试，匹配引擎根本没有提高；模式的剩余部分还要尝试界定符的右边。
(?! re)|	前向否定界定符。与肯定界定符相反；当所含表达式不能在字符串当前位置匹配时成功
(?> re)|	匹配的独立模式，省去回溯。
\A	|匹配字符串开始
\Z	|匹配字符串结束，如果是存在换行，只匹配到换行前的结束字符串。
\z	|匹配字符串结束
\G	|匹配最后匹配完成的位置。
\b	|匹配一个单词边界，也就是指单词和空格间的位置。例如， 'er\b' 可以匹配"never" 中的 'er'，但不能匹配 "verb" 中的 'er'。
\B	|匹配非单词边界。'er\B' 能匹配 "verb" 中的 'er'，但不能匹配 "never" 中的 'er'。
\n, \t, 等.	|匹配一个换行符。匹配一个制表符。等
\1...\9	|匹配第n个分组的内容。
\10	|匹配第n个分组的内容，如果它经匹配。否则指的是八进制字符码的表达式。


### 4.2 正则表达式实操学习

https://regexone.com/

#### 正则表达式练习
- 应用 notepad++ 或者editplus 软件 正则表达式处理文本
- 不足之处->无法处理大文件（内存限制）

##### 练习要求：
- 处理要求1：把序列转换成一行；
- 处理要求2：把一行的序列文件再转换回来；



处理文件内容示例：

```
>NP_651973.4 nbs, isoform E [Drosophila melanogaster]
MFVLTKDDEKFVLFPGKKVYTIGRLATDLIVAQDLSISRNHAQLLIQTEADGDDTLHIEDLGSRYGTFIF
PKNSQKPRKVPAKTSTPLPVGTRLRFGANMSIWQVTQLKLVTTVSALTRSEVQELTKMLEPMGGTVTSNW
TEECSHLTMNEVSVTVKLLHAMLENKPIVTFPYWRKMLQAAQSIHVKEGWPQPEDYQPTNIDVTWRPERT
RLFAGKTFVFMNRKHFDMYGSVVQKAGATCKDINSGVRKTFLTKSDVIVIQYVPSSQSQATESINSIQDR
YILEQNGRRIIQEYEIGMALIHCSITEFCNPTHKFISDSLPTTESVTSSMAFNSSIIVPNTERHSAQSNA
TPISELVVPESIECEMEQDASKPHSEDQASLRKRSHASTVDSSDEEKKSTLSKRAKSDIATKLTMKSKNA
ILLDSSLEEDVTPAPAPAPVQRVTRQSKAIAEEKSVHPPVPAASKHITRKTKQVFCVDSSDEENENARKP
KETPAPTIPSMAKKKTEAPVATRISPRLNGKSLATNITNQPADKHAVPAKRPVLSVASSDEEDEGDLFQF
RKSPQKPAETVVQPRIAGKGNAPARISVVDFLEKSQAQEPAPVPPQLESQSQTQPRKRLRLELLNESDSD
DCDNLFNFADSKKKRKTQEAQRNDDSTDGLFNFNSERPSDHDDEDSRLTEPFVPETESKKQSKYIVAPRR
DRPKKVDISGWLSCSRLNDNIKSEIDADSVKMETSIKADPDEEQWLAAMKDSIEVRMCNLNIVIRSQEEV
DASLEDSVNKHGGRKNFKKFVKTKNPHPQKRIVALKSLRLADGMVTCV

>AIQ85043.1 NBS-LRR disease resistance protein, partial [Musa ABB Group]
MGGVGKTTLAQQAYNPERVKDYFHHKVWVCVSDNFNVERLSKEIIESITENKCDLSNLDTLQVVVKKKLT
SKRFLLVLDDVWNEDSLKWERFCAPLRYGEPGSKILVTTRSKKIAEMVGNPFPLGGLDEASYWKLFKKCA
FGSEYAGE

>AIQ85044.1 NBS-LRR disease resistance protein, partial [Musa laterita]
GGGGKTSLAQQAYNHERVKDYFHHKVWVCVSDNFNVERLTKEIIESLTRNKWDLNNLDTLQVVVKEELTS
KRFLLVLDDVWNEDSLKWERFCAPLRYGEPGSKILVTTRSKKIAEMVGNPIPLGGLDEASYWELFKKCAF
GSEDAGE

>AIQ85045.1 NBS-LRR disease resistance protein, partial [Musa laterita]
MGGVGKTTLAQQAYNHERVQDYFQHEVWVCVSDNFNVERLTKEIIESITENKCDLSNLDTLQVVLKKNLT
SKRFLLVLDDVWNEDSLKWERFCAPLRYGEPGSKILVTTRSKNVFENGWNPIPLGGLDEASYWKLFKKCA
FGSEDAGEFPHLE

>AAM28915.1 NBS, partial [Pinus taeda]
TRFDWKEQLHRLQHVLPSETQEKLXFGYLNLNREERQMFLDSACFFIGQKRDTAIRIWEGSLWDGHSGFL
TLQHRCLLGVDDENNIEXHDHLRDFGRAACPNRFLPSWIPMDSLRVLQVSGSVLKTLWEDDSQPPLQLRE
LEINAPLSNIPGSIGRLKHLERFVVGKYLSGQVNLTELPVEFCHLQSLKALVLTECSKIKSLPEFGALLM
WLRHIDLSFCRNLERLPDSLHYLSHLRLINLSDCHDLVTLPDNIGRLRCLQHIDLQGCHNLERLPDSFGE
LTDLRHINLSGCHDLQRLPDSFGKLRYLQHIDLHGCHSLEGLPISFGDLMNLEYINLSNCHNLERLPESI
GNLSDLRHIDLSGCHNLERLPDNFRELEELRYLDVEGCSNLIIDRFEIIGISDNLPVAHQVNWNKY

```

### 4.3 Python中正则表达式re包


==新包学习思路总结：==

1. 学习里面的函数/方法
1. 类实例化后对象的学习(方法，属性)
1. 都有哪些类，类之间的关系

#### re包中常用方法



1. re.search(pattern, string, flags=0)
1. re.match(pattern, string, flags=0)
1. re.findall(pattern, string, flags=0)

参数说明
- pattern	匹配的正则表达式
- string	要匹配的字符串。
- flags	用于控制正则表达式的匹配方式，如：是否区分大小写，多行匹配等等。见：正则表达式修饰符 

### 注意区别：
**re.match与re.search的区别**
re.match只匹配字符串的开始，如果字符串开始不符合正则表达式，则匹配失败，函数返回None；而re.search匹配整个字符串，直到找到一个匹配。


**re.match和re.search 与findall的区别：**
前面两个找到一个就结束，findall会找到所有。
返回对象不同，findall返回列表，


#### 正则表达式修饰符 - 可选标志
正则表达式可以包含一些可选标志修饰符来控制匹配的模式。修饰符被指定为一个可选的标志。多个标志可以通过按位 OR(|) 它们来指定。如 re.I | re.M 被设置成 I 和 M 标志：

修饰符|	描述
--|--
re.I|	忽略大小写
re.M|	多行匹配，影响 ^ 和 $
re.U|	根据Unicode字符集解析字符。这个标志影响 \w, \W, \b, \B.
re.L|	做本地化识别（locale-aware）匹配
re.S|	使 . 匹配包括换行在内的所有字符
re.X|	该标志通过给予你更灵活的格式以便你将正则表达式写得更易于理解。

代码示例：

```
import re
 
line = "Cats are smarter than dogs"
 
matchObj = re.match( r'dogs', line, re.I)
if matchObj:
   print("match --> matchObj.group() : ", matchObj.group())
else:
   print ("No match!!")


line = "Cats are smarter than Dogs"
matchObj = re.search( r'dogs', line, re.I)
if matchObj:
   print ("search --> matchObj.group() : ", matchObj.group())
else:
   print ("No match!!")
```



#### 匹配对象的方法  得到详细的搜索匹配结果

- group() 返回被 RE 匹配的字符串。
- start() 返回匹配开始的位置
- end() 返回匹配结束的位置
- span() 返回一个元组包含匹配 (开始,结束) 的位置

更多方法及帮助：https://docs.python.org/zh-cn/3.7/library/re.html#match-objects

```
import re


line = "Cats are smarter than dogs";
 
searchObj = re.search( r'(.*) are (.*?) .*', line, re.M|re.I)
 
if searchObj:
   print "searchObj.group() : ", searchObj.group()
   print "searchObj.group(1) : ", searchObj.group(1)
   print "searchObj.group(2) : ", searchObj.group(2)
else:
   print "Nothing found!!"
   
   
print(re.search('www', 'www.omicsclass.com').span())  # 在起始位置匹配
print(re.search('com', 'www.omicsclass.com').span())         # 不在起始位置匹配   
   
```


#### 搜索替换
Python 的 re 模块提供了re.sub用于替换字符串中的匹配项。支持正则表达式，比字符串自带的replace方法功能更强大。

语法：

>re.sub(pattern, repl, string, count=0, flags=0)

- pattern : 正则中的模式字符串。 (搜索的内容)
- repl : 替换的字符串，也可为一个函数。 (替换的内容)
- string : 要被查找替换的原始字符串。 (在哪里搜索替换)
- count : 模式匹配后替换的最大次数，默认 0 表示替换所有的匹配。
- flags	用于控制正则表达式的匹配方式，如：是否区分大小写，多行匹配等等。见：正则表达式修饰符 


```
import re
 
tel = "010-8054-3251,www.omicsclasss.com"

# 删除字符串中的 数字 
num = re.sub(r"^\d+-\d+-\d+,", "", tel)
print("Website is:", num)
 
# 删除数字
num = re.sub(r'\d+', "", tel)
print("telphone num : ", num)



### 高级用法 （选修）
#repl 参数是一个函数
#以下实例中将字符串中的匹配的数字乘以 2：

# 将匹配的数字乘以 2
def double(matched):
    value = int(matched.group('value'))
    return str(value * 2)
 
s = 'A23G4HFD567'
print(re.sub('(?P<value>\d+)', double, s))
```
#### re.split 方法支持正则表达式：
split 方法按照能够匹配的子串将字符串分割后返回列表，它的使用形式如下：

>re.split(pattern, string[, maxsplit=0, flags=0])
- pattern	匹配的正则表达式
- string	要分割的字符串。
- maxsplit	分隔次数，maxsplit=1 分隔一次，默认为 0，不限制次数。
- flags	标志位，用于控制正则表达式的匹配方式，如：是否区分大小写，多行匹配等等。


```
import re
re.split(r'[,;=]', 'omicsclass, omicsclass;omicsclass.')
re.split(r'\s+', ' omicsclass, omicsclass, omicsclass.') 
re.split(r',', ' omicsclass, omicsclass, omicsclass.', 1) 

 
re.split(r'\t', 'hello world')   # 对于一个找不到匹配的字符串而言，split 不会对其作出分割
['hello world']

```


#### re.compile 方法
compile 方法用于编译正则表达式，提前编译后的正则表达式可以加快正则表示的匹配速度；编译后生成一个正则表达式（ Pattern ）对象。

语法格式为：

>re.compile(pattern[, flags])

- pattern : 一个字符串形式的正则表达式
- flags	用于控制正则表达式的匹配方式，如：是否区分大小写，多行匹配等等。见：正则表达式修饰符 


#####  Pattern正则表达式对象 （正则对象）

常见方法：
- Pattern.findall(string[, pos[, endpos]])
- Pattern.search(string[, pos[, endpos]])
- Pattern.match(string[, pos[, endpos]])

参数：
- string : 待匹配的字符串。
- pos : 可选参数，指定字符串的起始位置，默认为 0。
- endpos : 可选参数，指定字符串的结束位置，默认为字符串的长度。

注意： match 和 search 是匹配一次 findall 匹配所有。

更多方法：https://docs.python.org/zh-cn/3.7/library/re.html#regular-expression-objects

```
import re
 
p = re.compile(r'\d+')   # 查找数字
result1 = p.findall('omicsclass 123 google 456')
result2 = p.findall('runomicsclass123google456', 0, 10)
 
print(result1)
print(result2)



result3=p.search("runomicsclass123google456")
result4=p.match("runomicsclass123google456")

print(result3)
print(result4)

#或者也可以这样用
result5=re.search(p,"runomicsclass123google456")
result6=re.match(p,"runomicsclass123google456")

print(result5)
print(result6)
```



#### 正则表达式处理文件主要用处

1. split分隔文件
2. 查找是否含有某某字符，配合IF判断
3. 在字符串中捕获信息


#### 练习

从GFF文件中提取基因ID及位置信息存成表格(注意用正则表达式的捕获功能完成)

GFF 文件（人1号染色体上基因注释信息）地址：ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.chromosome.1.gff3.gz

答案：

```
import re

fr=open("D:\\python_script\\Homo_sapiens.GRCh38.98.chromosome.1.gff3\\Homo_sapiens.GRCh38.98.chromosome.1.gff3","r")

fw=open("D:\\python_script\\chr1.txt","w")

for line in fr:
    if re.match("#",line):
        continue
        
    tmp=re.split("\t",line)
    if tmp[2] == "gene":
        mobj=re.search("ID=gene:([^;]+)",tmp[8]) #捕获基因ID
        if mobj:
            fw.write("\t".join([tmp[0],tmp[3],tmp[4],tmp[6],mobj.group(1)])+"\n")
            
fr.close()
fw.close()
        
    
    
```


---

## 5 Python biopython包处理生物数据


### 5.1 BioPython包



专门用于处理生物数据的包，主要功能如下：
#### （1）biopython处理生物数据
- Blast output – both from standalone and WWW Blast
- Clustalw
- FASTA/FASTQ
- GenBank
- PubMed and Medline
- ExPASy files, like Enzyme and Prosite
- SCOP, including ‘dom’ and ‘lin’ files
- UniGene
- SwissProt

#### （2）在线链接生物数据库，实时处理数据

- NCBI – Blast, Entrez and PubMed services
- ExPASy – Swiss-Prot and Prosite entries, as well as Prosite searches

#### （3）与生物信息常用软件交户，实现批量处理分析数据

- Standalone Blast from NCBI
- Clustalw alignment program
- EMBOSS command line tools


#### 获取帮助：
官方网站：https://biopython.org/

帮助文档：http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

中文帮助文档（Last Update – 22 March 2013 (Biopython 1.61+)）：https://biopython-cn.readthedocs.io/zh_CN/latest/

#### 安装


```
pip install biopython -i  https://pypi.tuna.tsinghua.edu.cn/simple #安装
pip install --upgrade biopython  #更新
pip uninstall biopython  #卸载

#导入包：
import Bio
```
### 5.2 利用biopython处理序列(fasta，fastq等)
#### Seq序列对象与SeqRecord注释对象学习
**Seq对象中文帮助：** https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr03.html#chapter-bio-seq

**SeqRecord注释对象中文帮助：** https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr04.html#chapter-seqrecord


#### 5.2.1 Seq序列对象
```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

#Seq 对象创建
dna_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
protein_seq = Seq("EVRNAK", IUPAC.protein)

#序列对象继承了字符串对象一些方法
len(dna_seq)
dna_seq.count("A")
GC(dna_seq)  #计算GC含量
dna_seq[4:12] #切取序列
str(dna_seq)  #将序列对象转换成字符串
my_seq = Seq("ACGT", IUPAC.unambiguous_dna)
my_seq + dna_seq  #连接或添加序列
my_seq.upper()  
my_seq.lower()

#核苷酸序列反向互补序列            
my_seq.reverse_complement()  

#转录
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
messenger_rna = coding_dna.transcribe()
messenger_rna.back_transcribe()  #从mRNA逆向转录为DNA编码链的方法

#翻译  https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr03.html#sec-translation
coding_dna.translate()

```


Bio.Alphabet.IUPAC 提供了蛋白质、DNA和RNA的基本定义：

- IUPAC.protein  蛋白质
- IUPAC.unambiguous_dna  DNA

更多编码见：https://www.omicsclass.com/article/409




#### 5.2.2 SeqRecord序列注释对象
SeqRecord 类非常简单,包括下列属性:


**.seq**
– 序列自身（即 Seq 对象）。

**.id**
– 序列主ID（-字符串类型）。通常类同于accession number。

**.description**
– 序列描述（-字符串类型）。

**.name**
– 序列名/id （-字符串类型）。 可以是accession number, 也可是clone名（类似GenBank record中的LOCUS id）。

**.letter_annotations**
– 对照序列的每个字母逐字注释（per-letter-annotations），以信息名为键（keys），信息内容为值（value）所构成的字典。值与序列等长，用Python列表、元组或字符串表示。.letter_annotations可用于质量分数(如第 18.1.6 节) 或二级结构信息 (如 Stockholm/PFAM 比对文件)等数据的存储。

**.annotations**
– 用于储存附加信息的字典。信息名为键（keys），信息内容为值（value）。用于保存序列的零散信息（如unstructured information）。

**.features**
– SeqFeature 对象列表，储存序列的结构化信息（structured information），如：基因位置, 蛋白结构域。

**.dbxrefs**
– 储存数据库交叉引用信息（cross-references）的字符串列表。

#### SeqRecord对象创建
```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

simple_seq = Seq("CCCTTCTTGTCTTCAGCGTTTCTCC", IUPAC.unambiguous_dna)
simple_seq_r = SeqRecord(simple_seq, id="AC12345",description="just a test sequence")

#属性后添加和修改
simple_seq_r.description = "just a test sequence"

simple_seq_r.seq
simple_seq_r.id
simple_seq_r.description

#如果是fastq文件，可以添加序列质量值注释信息

simple_seq_r.letter_annotations["phred_quality"] = [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26,
26, 26, 26, 23, 23]
print(simple_seq_r.letter_annotations)
print(simple_seq_r.letter_annotations["phred_quality"])


#序列截取：

simple_seq_r[0:2]

```

#### SeqRecord对象从文件中获得
fasta/fastq文件读写，得到SeqRecord对象；

```
#fasta文件读取得到SeqRecord对象
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

os.chdir("D:\\python_script")
output_handle = open("out.fa", "w")
for rec in SeqIO.parse("test.fa", "fasta"):
    seq=rec.seq
    seq_r = SeqRecord(seq[0:10],id=rec.id,description=rec.description)
    SeqIO.write(seq_r, output_handle, "fasta")

output_handle.close()
```
#### 5.3 练习
##### 任务：利用biopython处理fasta/fastq序列

1. 任务: 提取指定ID的fasta序列
1. 任务: 截取fasta序列中指定位置的序列
1. 任务: fastq文件转换成fasta文件
1. 任务: fastq文件去掉前5个碱基

fa文件：
```
>AT1G66550.1 cds chromosome:TAIR10:1:24828537:24829589:1 gene:AT1G66550 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:WRKY67 description:Probable WRKY transcription factor 67 [Source:UniProtKB/Swiss-Prot;Acc:Q93WV7]
ATGGTTTCCAACATTGATCACAAGGCTATGGAAGCACTCCTCCGTGGCCAAGGATGCGCT
AACAACCTCAAGATTCTCCTTGAAAACGGCGAAATAAGCTCAGTTTCAACAGAACCACTC
ATCCACACCATTCTCGATTCTTTCTCACTTGCTCTGTCTTTTATGGATTCTCCTAATCAT
CCACCATACCATGAATCCTCTTCTCATAACATGGCAAGTCATATGTCCCGGAGATCATCT
AAGCAAGTACAACATCGCAGAAAACTTTGTGTAGCAGAAGGTTTAGTGAATTACAATCAC
GATTCCCGGACTATGTGCCCCAATGATGGCTTCACCTGGAGGAAATATGGACAAAAAACC
ATTAAAGCCTCAGCGCACAAAAGGTGTTACTATCGGTGTACCTATGCAAAAGACCAAAAC
TGCAATGCTACAAAGCGGGTGCAGAAGATCAAAGACAACCCTCCAGTGTACAGAACCACT
TACTTGGGAAAACATGTGTGTAAAGCTTTTGCAGTTCATGATGATACATATAGTTCCACG
ATGATTCGATTCGACCAAGTTGTTCCTGAACCGATTATGCCGCAGCTCACAACAATTGAC
CACCAAGTAATTACCGTGGAGGAAAACTCCGCAGAACATATCATGAACCAAGAATGTGAT
ATTAATGATTATTTGGTGGATGATGACCCATTTTGGGCTAGTCAATTTCCCCCGTTTCCA
TCGAGTGACACAATGTTCTTGGAAAACATTTCTGCTTTTGATTAG
>AT4G39410.1 cds chromosome:TAIR10:4:18332606:18334893:-1 gene:AT4G39410 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:WRKY13 description:WRKY13 [Source:UniProtKB/TrEMBL;Acc:A0A178UV80]
ATGGGTGCGATAAACCAAGGAATAAGCTTGTTTGATGAATCACAAACCGTCATAAACCCT
ATTAATACCAACCATCTAGGTTTCTTCTTCTCTTTCCCTAGTCACAGCACCTTATCTTCA
TCATCTTCGTCGTCTTCGTCTTCTCCTTCTTCTCTTGTGTCTCCATTTCTTGGTCATAAC
TCCCTAAACTCCTTCCTTCATAATAACCCGTCTTCATTCATAAGTCATCCTCAAGATTCC
ATCAATCTCATGACCAATCTCCCCGAAACCCTAATCTCGTCTTTGTCCTCATCAAAGCAA
AGGGACGATCATGATGGTTTTCTTAATCTCGATCATCATCGTCTTACCGGTAGTATTTCA
TCCCAAAGACCCCTGTCAAATCCATGGGCATGGAGTTGTCAAGCGGGATACGGAAGCAGC
CAGAAAAACAACCATGGAAGCGAGATTGATGTTGATGATAATGATGATGAGGTTGGCGAT
GGTGGTGGCATTAATGATGATGATAATGGTCGTCATCATCATCATGATACTCCCAGTCGT
CATGATAAACATAACACAGCGTCATTAGGCGTAGTTTCTTCTCTGAAGATGAAGAAGCTT
AAGACAAGAAGAAAAGTGAGGGAACCTCGGTTTTGCTTTAAGACACTTAGCGAGGTTGAT
GTCTTAGATGATGGATATAGATGGAGAAAGTATGGCCAGAAAGTTGTCAAAAACACCCAA
CATCCCAGGAGCTATTACAGATGCACACAAGACAAGTGTAGAGTGAAGAAGAGAGTGGAG
AGATTAGCAGATGACCCAAGAATGGTAATCACTACTTACGAAGGAAGACACCTTCACTCT
CCTTCTAATCATCTCGACGACGACTCTCTCTCCACCTCTCACCTGCACCCTCCTCTCTCC
AACTTCTTCTGGTGA
>AT4G18170.1 cds chromosome:TAIR10:4:10061214:10062893:1 gene:AT4G18170 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:WRKY28 description:WRKY28 [Source:UniProtKB/TrEMBL;Acc:A0A178V3M3]
ATGTCTAATGAAACCAGAGATCTCTACAACTACCAATACCCTTCATCGTTTTCGTTGCAC
GAAATGATGAATCTGCCTACTTCAAATCCATCTTCTTATGGAAACCTCCCATCACAAAAC
GGTTTTAATCCATCTACTTATTCCTTCACCGATTGTCTCCAAAGTTCTCCAGCAGCGTAT
GAATCTCTACTTCAGAAAACTTTTGGTCTTTCTCCCTCTTCCTCAGAGGTTTTCAATTCT
TCGATCGATCAAGAACCGAACCGTGATGTTACTAATGACGTAATCAATGGTGGTGCATGC
AACGAGACTGAAACTAGGGTTTCTCCTTCTAATTCTTCCTCTAGTGAGGCTGATCACCCC
GGTGAAGATTCCGGTAAGAGCCGGAGGAAACGAGAGTTAGTCGGTGAAGAAGATCAAATT
TCCAAAAAAGTTGGGAAAACGAAAAAGACTGAGGTGAAGAAACAAAGAGAGCCACGAGTC
TCGTTTATGACTAAAAGTGAAGTTGATCATCTTGAAGATGGTTATAGATGGAGAAAATAC
GGCCAAAAGGCTGTAAAAAATAGCCCTTATCCAAGGAGTTACTATAGATGTACAACACAA
AAGTGCAACGTGAAGAAACGAGTGGAGAGATCGTTCCAAGATCCAACGGTTGTGATTACA
ACTTACGAGGGTCAACACAACCACCCGATTCCGACTAATCTTCGAGGAAGTTCTGCCGCG
GCTGCTATGTTCTCCGCAGACCTCATGACTCCAAGAAGCTTTGCACATGATATGTTTAGG
ACGGCAGCTTATACTAACGGCGGTTCTGTGGCGGCGGCTTTGGATTATGGATATGGACAA
AGTGGTTATGGTAGTGTGAATTCAAACCCTAGTTCTCACCAAGTGTATCATCAAGGGGGT
GAGTATGAGCTCTTGAGGGAGATTTTTCCTTCAATTTTCTTTAAGCAAGAGCCTTGA
>AT5G45050.1 cds chromosome:TAIR10:5:18176914:18181973:-1 gene:AT5G45050 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:RRS1B description:Probable WRKY transcription factor 16 [Source:UniProtKB/Swiss-Prot;Acc:Q9FL92]
ATGACCGAGAGTGAGCAAATCGTCTACATCAGCTGCATAGAGGAGGTACGATACTCCTTC
GTCAGCCACCTCTCCAAAGCTCTCCAGCGAAAAGGTGTAAACGATGTCTTCATCGATAGC
GATGATTCGCTTTCCAACGAGTCTCAATCAATGGTCGAGAGAGCTAGGGTTTCTGTTATG
ATTTTACCAGGAAACCGTACGGTATCTCTTGACAAGCTCGTGAAGGTTCTCGATTGCCAG
AAGAACAAAGATCAAGTGGTGGTTCCGGTGTTGTACGGTGTCAGATCATCAGAGACCGAA
TGGCTTAGCGCGCTGGATTCGAAAGGATTCTCATCAGTACACCATTCCAGGAAAGAATGT
AGTGACTCCCAGCTTGTAAAAGAGACTGTTAGAGATGTGTATGAGAAGCTCTTTTATATG
GAACGAATTGGAATATATTCGAAGCTGCTGGAGATTGAGAAAATGATTAACAAGCAACCG
TTGGACATCCGTTGTGTTGGAATTTGGGGTATGCCTGGCATAGGCAAGACTACACTTGCT
AAAGCAGTCTTTGACCAAATGTCTGGTGAGTTTGATGCTCATTGCTTTATTGAAGACTAC
ACCAAAGCTATTCAAGAGAAGGGTGTTTATTGTTTGCTGGAGGAACAGTTTTTGAAAGAA
AATGCTGGTGCTAGTGGTACCGTTACGAAATTGAGCTTGCTTAGGGATAGATTAAACAAT
AAGAGGGTTCTTGTTGTTCTTGATGATGTCCGCAGTCCTCTGGTTGTGGAGTCTTTTCTT
GGAGGGTTTGACTGGTTTGGTCCCAAAAGTCTAATCATCATAACCTCCAAAGATAAATCG
GTGTTTCGCCTTTGTCGAGTCAATCAAATATACGAGGTTCAGGGTTTAAATGAGAAAGAG
GCTCTTCAACTCTTCTCTTTGTGTGCGTCTATAGACGATATGGCAGAGCAGAATCTCCAC
GAGGTGTCAATGAAAGTTATTAAATATGCTAATGGCCATCCATTAGCTCTCAATCTCTAT
GGCAGAGAACTGATGGGGAAGAAAAGACCACCAGAAATGGAGATAGCATTCCTCAAACTC
AAGGAATGTCCTCCAGCTATTTTTGTTGATGCAATCAAGAGCTCGTATGACACACTCAAT
GACAGGGAAAAAAACATTTTTTTGGACATAGCTTGTTTCTTCCAGGGAGAAAATGTTGAC
TACGTGATGCAACTGCTTGAGGGTTGTGGTTTCTTTCCACATGTTGGAATTGATGTTCTT
GTGGAGAAGAGTCTGGTGACTATTTCAGAAAACCGAGTGCGGATGCATAACTTGATCCAA
GATGTTGGCCGACAAATAATAAATAGAGAAACAAGACAGACTAAGAGGCGCAGCAGACTG
TGGGAACCTTGCAGCATCAAATATTTATTAGAAGATAAGGAACAAAACGAAAATGAAGAA
CAAAAAACAACTTTTGAACGTGCTCAGGTCCCTGAAGAGATCGAAGGCATGTTTCTGGAC
ACATCAAACTTAAGTTTTGATATTAAGCATGTTGCCTTTGATAATATGTTGAACCTTAGA
TTGTTCAAGATTTACAGTTCCAATCCTGAAGTCCATCATGTAAACAATTTCCTCAAAGGC
TCTCTCAGTTCTCTTCCTAATGTGCTAAGACTCCTGCATTGGGAGAACTATCCTCTGCAG
TTTCTGCCTCAAAATTTTGATCCTATACACCTTGTTGAAATCAACATGCCGTACAGCCAA
CTTAAGAAACTTTGGGGTGGAACCAAGGACCTGGAGATGTTGAAGACAATCAGGCTTTGT
CATTCCCAACAACTAGTTGATATTGACGATCTTTTAAAAGCTCAAAATCTTGAGGTAGTT
GATCTCCAAGGCTGTACAAGACTGCAGAGTTTCCCAGCCACCGGTCAATTGCTACATTTA
CGAGTTGTAAATCTCTCAGGTTGCACAGAGATCAAAAGTTTCCCAGAAATTCCCCCAAAT
ATTGAGACACTGAATCTACAGGGGACTGGTATAATAGAATTACCACTTTCCATTGTTAAG
CCAAACTACAGAGAGCTTTTGAATCTTCTAGCTGAAATCCCGGGTCTTTCAGGTGTCTCA
AACCTTGAGCAAAGTGATCTCAAACCTTTAACAAGCCTGATGAAAATTAGCACATCTTAC
CAAAATCCTGGCAAGCTTAGTTGCTTGGAGCTGAATGATTGTTCTCGTTTGCGAAGTCTG
CCAAACATGGTTAATTTAGAACTTCTCAAAGCCCTTGATCTTTCTGGTTGCTCAGAGCTC
GAGACTATCCAGGGTTTCCCACGGAACCTGAAAGAGTTATATCTTGTTGGCACTGCAGTA
AGACAAGTGCCACAACTTCCTCAAAGTCTAGAATTCTTTAATGCCCATGGTTGTGTCTCT
CTCAAATCAATTCGTTTGGACTTCAAGAAGCTTCCTGTGCATTACACATTTAGTAATTGT
TTCGATCTATCTCCACAAGTGGTCAACGATTTTTTAGTGCAGGCGATGGCTAATGTGATT
GCAAAACACATACCAAGAGAGCGTCATGTCACAGGCTTTTCTCAAAAGACTGTGCAGCGT
TCGAGTCGTGACAGTCAGCAGGAACTCAACAAAACTTTGGCTTTCAGCTTCTGTGCGCCC
TCACATGCGAATCAAAATTCCAAACTTGATCTGCAACCAGGATCTTCTTCAATGACACGA
CTAGATCCTTCTTGGAGGAACACACTTGTGGGCTTTGCTATGCTGGTGCAAGTCGCATTT
TCCGAGGGTTACTGTGATGATACTGATTTTGGCATTAGTTGTGTTTGCAAATGGAAAAAC
AAGGAAGGCCACTCTCATAGGAGAGAAATAAATTTGCATTGTTGGGCTTTAGGGAAAGCT
GTTGAAAGGGATCATACGTTTGTCTTCTTTGATGTCAACATGCGTCCAGATACCGATGAA
GGAAATGACCCCGATATCTGGGCTGATTTAGTTGTTTTTGAGTTCTTTCCTGTCAATAAA
CAGAGAAAGCCTCTAAATGATAGTTGCACAGTGACAAGATGTGGAGTCCGTTTAATAACT
GCTGTAAACTGCAATACAAGTATCGAGAATATATCACCAGTTTTGTCCTTGGATCCGATG
GAGGTTTCTGGTAATGAAGATGAAGAAGTATTGAGAGTCAGATATGCTGGTTTACAGGAG
ATATATAAAGCTTTGTTTCTTTACATAGCGGGTTTGTTCAATGACGAGGATGTTGGTTTG
GTAGCACCACTTATTGCTAACATTATTGACATGGACGTTAGTTATGGGCTCAAGGTCTTA
GCCTATAGGTCTCTCATACGTGTATCTTCCAATGGGGAAATAGTGATGCACTATTTGCTA
CGACAAATGGGTAAAGAAATCCTCCATACAGAATCAAAGAAGACTGACAAATTAGTCGAC
AATATTCAGAGTTCCATGATCGCAACAAAGGAAATCGAGATCACTCGTTCAAAGAGTCGC
CGAAAGAACAACAAGGAAAAGAGAGTGGTTTGCGTAGTGGATCGAGGCAGCCGGTCCAGT
GACCTATGGGTTTGGCGAAAGTATGGTCAAAAACCCATCAAAAGTTCTCCTTATCCAAGG
AGTTACTATAGATGTGCCAGCTCGAAAGGTTGTTTTGCTAGGAAACAAGTCGAACGTAGC
CGCACTGATCCAAATGTTTCAGTAATTACTTACATCTCTGAGCATAACCATCCATTCCCC
ACTCTACGCAATACTCTTGCCGGCTCCACTCGTTCCTCTTCCTCCAAATGCTCAGATGTA
ACTACTTCTGCCTCATCGACAGTCTCCCAAGACAAAGAAGGACCGGATAAATCCCATTTG
CCTTCCTCCCCTGCTTCTCCTCCTTATGCGGCCATGGTGGTTAAGGAGGAGGACATGGAG
CAATGGGACAATATGGAGTTCGATGTTGACGTTGAAGAAGATACTTTCATACCCGAATTA
TTTCCAGAGGATACCTTCGCTGATATGGACAAGCTTGAGGAAAATTCTCAGACTATGTTT
CTCTCTCGCAGAAGCAGCGGAGGCAACATGGAAGCCCAAGGGAAGAACTCTAGTGATGAT
AGGGAGGTCAATTTACCTAGTAAAATTCTGAATAGATAG

```

fastq文件 ：
详细介绍可以观看：https://www.omicsclass.com/course/28

```
@A00808:122:HM2JCDSXX:1:1101:3025:1000 1:N:0:TCGGAAGT+ACAAGGCT
AATTATTCCAAGTGATCATTTAATTTAACAGTACGTTATTACAGTTTTTGACAATACACCAGGAGGGGCAGAAGCAGCACTCCATTTATGCGCCGAATTCTCGTACAGACACACACACACACACAGATTCACACAGTCAATTCAGTCCTT
+
FFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:
@A00808:122:HM2JCDSXX:1:1101:4399:1000 1:N:0:TCGGAAGT+ACAAGGCT
CTTTGCTATAAACTAAGTCACCTTCTACAGCCTGCGAAATGCCATATTTTTCAACTCTGGCACTGGCAGCATGGTTCCAGAGGTAACTTTGGTAACTATGAACGTACATCATTCGTAATGTTCTTGGTATACTCTTTAATGCTTGCAGAT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00808:122:HM2JCDSXX:1:1101:8522:1000 1:N:0:TCGGAAGT+ACAAGGCA
NACTTTTCAGATTTCAACTTGAGAGATGGACCATCAAGAGCCTTTTTCATGTCACGTTTTTGAGCCATACAATGATACAGACGAAAGAAAACAAGCCAAACAACTCATCCAACCGACCGCGAACAGTAGACACAATTACACAAACACATG
+
#FFF,:FFF:,:FFF:F:FF:FFFFFF:FFF,F:FFF,F,:FF,:FFFF:FFFFF:FF,FF::,,FFFFFFFF:F:FFF,,,FFFFFF,FFF,:,:F::FFFFFFFFFFFFFF,FFFFFF:FFF,,F:,,F:F,:,FFF:,,FF:,FFFF
@A00808:122:HM2JCDSXX:1:1101:13205:1000 1:N:0:TCGGAAGT+ACAAGGCT
NTTCGACGAGGAGGCCGCCTGCGCGCGCGACGCCGCCGGGGAGGCCCTCGCGGCCTTCGAGTCGCTGCTCGCGCGCCTCCCCCCGCCCGACGCCGACTCGCGCCGCCGATCCATGGGGCTCAAGTTGGAGCAGCACAAAGCCGAGCTCAA
+
#FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF,F:FFFF,FFFFF,F,FFFFFF:FF:FF:FF:,:,,F,FFFFF,F,:F:FFF:,F,FFFFF,F,:,,,FFF::,:FFFF,FF:,:,:,,FFF,F,FF,,FFFFF,FFFF,FFFF
@A00808:122:HM2JCDSXX:1:1101:3902:1016 1:N:0:TCGGAAGT+ACAAGGCT
GTCGTCGTCACCGAATGTTGTCTTCTTTCCTGCACTAGGTGTACCAGCGCTCTGCTTGAACGGTGTACCCCGTCCACCACGGCCCCTGTCACCCCTACCAAAGCCTCTACCACGTCCACGGTCACCGCGCCCACGGCTACCATCACCACG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF
@A00808:122:HM2JCDSXX:1:1101:20518:1016 1:N:0:TCGGAAGT+ACAAGGCT
NGCTGGAGAACTAAGTGTCGATCTTCACCTTCACCAACAGATTCATGAACAAACCTAGCAGCGACACAAAAGCAAGACATTAGTAGTTGCAAATATTACTCTATCCAGTTAAGTTAGTATTTGAGAAGAGGGTAGGTGAAGCAGGTCGAA
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00808:122:HM2JCDSXX:1:1101:23158:1016 1:N:0:TCGGAAGT+ACAAGGCT
NTCCCTGGACTCTGGGGCACTTTAGGGCTAGGAACACTTGTCAAGCGACTCTGAGGCACAGTATCAGACATGGACTTCGAGAAGGCTTGTGCAAGCTCCAATGCAGATGCTCCTTTCCCAAAACGAGGTGCAACAACCTCAGGCTTTGGT
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF
@A00808:122:HM2JCDSXX:1:1101:1597:1031 1:N:0:TCGGAAGT+ACAAGGCT
CGCCGCCGCACTCTGCAACCGCGCCGCGGGGAGGGAGGGAAGGACGAAGGAGGAAGGAGATGGAGCGGGTCGGCGGCGGGGAGAAGCAGCTGGAGGACTGCACCGTGTCCAATGCTCTCGGCACCTGGTTCTTCTCAGTTGCTGGTGCTC
+
F:FFFFFFFFFFFFFFFFFFFFFFFFFF,FF:FFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF
@A00808:122:HM2JCDSXX:1:1101:3351:1031 1:N:0:TCGGAAGT+ACAAGGCT
CCCCCTTGAACTATCCGTAAATACCACCATTCCACCAAACTCCCCTCTGTCTTCTACGACAAATAATCTTTTCAACCCAGTAACTTATAATTAATTACAGAGAGATTAGCAAACATATAATGGTAAGAATCAGTTGTCGTAGGTGCGGCA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF
@A00808:122:HM2JCDSXX:1:1101:5466:1031 1:N:0:TCGGAAGT+ACAAGGCT
CTCGCAGTGCTTGGTGTACTTGGCCGGGTCCAGCAGGGCGTCCAGCTCCGACGGGCTGCTCGGCTTCACCTTGATCATCCACCCGTCCTCGTACGGGCTTGAGTTAATCAGGCCGGGTGTCTCAGAGAGCTTGTCGTTAACCTCGACGAC
+
:FF:FF:F,FFFFF,FFFF:::FFF,F::F:F:,F,,,FFF:F::FFFF:FFF:F,:::FF:F,:FFFFFFFFFF,FFFFFFF::FFFF:FFFFFFFF:,FFFFFFFF,F,FF,F,FFF,:F,FFF,FFFF,FFFFF,F:,FFFFFFFFF
@A00808:122:HM2JCDSXX:1:1101:6714:1031 1:N:0:TCGGAAGT+ACAAGGCT
CGACGAGATCGGCCGCTCCGAATCCGATCCGCGCACGTGGTCTCGCCGGCGGGGTCAGGTTGGTTCGCGGGAGCCGCCGCCGCCGCCGACGACGACGACGATGGCCATGGAGACGCCGCCGCCGTTCCAGGAGTCCGCCCACTGCGACGT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

```

==任务1 答案：==
```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

os.chdir("D:\\python_script")

idlist={}  #字典用于存储ID列表
f = open("idlist.txt", "r")  #打开ID列表文件
for line in f:
    line=line.strip()
    idlist[line]=1
f.close()

f_out = open("get.fa", "w")
for rec in SeqIO.parse("test.fa", "fasta"):
    if rec.id in idlist:  #判断ID是否存在与ID列表字典中
        SeqIO.write(rec,f_out,"fasta")   #如果存在写出该序列
        
f_out.close()
```
==任务2 答案：==

```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import re
os.chdir("D:\\python_script")

idlist={}  #字典用于存储ID列表以及位置
f = open("id_pos.txt", "r")  #打开ID列表文件
for line in f:
    line=line.strip()
    tmp=re.split(r"\t",line)
    idlist[tmp[0]]=[int(tmp[1]),int(tmp[2])]
f.close()

f_out = open("get_pos.fa", "w")
for rec in SeqIO.parse("test.fa", "fasta"):
    if rec.id in idlist:  #判断ID是否存在与ID列表字典中
        start=idlist[rec.id][0]  #取得对应ID要截取的起始位置
        end=idlist[rec.id][1]   #取得对应ID要截取的结束位置
        rec_new=SeqRecord(rec.seq[start-1:end],id=rec.id,description=rec.description) #注意生物数据一般是从1作为索引，编程语言一般是从0开始所以要减一
        SeqIO.write(rec_new,f_out,"fasta")   #写出截取好的序列
        
f_out.close()
```

==任务3 答案：==
```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

os.chdir("D:\\python_script")
output_handle = open("fq2fa.fa", "w")
for rec in SeqIO.parse("test.fq", "fastq"):
    SeqIO.write(rec, output_handle, "fasta")

output_handle.close()
```


==任务4 答案：==

```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

os.chdir("D:\\python_script")
output_handle = open("trimed.fq", "w")
for rec in SeqIO.parse("test.fq", "fastq"):
    rec_new=rec[5:]
    SeqIO.write(rec_new, output_handle, "fastq")
output_handle.close()
```

## 6  Python 中数据统计分析包numpy 与 pandas

### 6.1 数据统计分析包介绍

- **NumPy：N维数组矩阵容器；** 基础的数学计算模块。

- **Pandas：表格容器；** 提供了一套名为DataFrame的数据结构，适合统计分析表格类数据。

- **SciPy：科学计算函数库；** 基于Numpy，提供方法(函数库)直接计算结果；封装了一些高阶抽象的物理模型。比方说做个傅立叶变换，做个滤波器等等。

==非数学研究，建议直接入手Numpy，pandas。目的，是方便我们处理表格类的生物数据。==


### 6.2 Numpy包学习

Numpy：
用来存储和处理大型矩阵，本身是由C语言开发，和列表（list）的使用很像，但是计算速度要比列表快很多。

- 官方网站：https://numpy.org/
- 中文文档：http://liao.cpython.org/parttwo/

#### 6.2.1 安装

```
# 使用 pip
pip install numpy -i  https://pypi.tuna.tsinghua.edu.cn/simple

#导入：
import numpy as np
```

#### 6.2.2 矩阵/数组创建函数


函数|描述
--|--
np.array()|创建数组
np.zeros()|创建数据全为0
np.ones()|创建数据全为1
np.empty()|创建数据接近0
np.arange()|按指定范围创建数据
np.linspace()|创建线段

```
import numpy as np #为了方便使用numpy 采用np简写
#列表转化为矩阵
a = np.array([2,23,4])
a = np.array([[11, 12, 13, 14, 15],
              [16, 17, 18, 19, 20],
              [21, 22, 23, 24, 25],
              [26, 27, 28 ,29, 30],
              [31, 32, 33, 34, 35]])
print(a)

#指定数据类型
a = np.array([2,23,4],dtype=np.int)
print(a.dtype)

a = np.arange(10,20,2) # 10-19 的数据，2步长
a = np.linspace(1,10,20)    # 开始端1，结束端10，且分割成20个数据，生成线段

a = np.zeros((3,4)) # 数据全为0，3行4列
a = np.ones((3,4),dtype = np.int)   # 数据为1，3行4列
a = np.empty((3,4)) # 数据为empty，3行4列
```

##### 数据类型

NumPy支持的数值类型

符号|	含义
--|--
np.bool|	True和Flase
np.int|	支持int的32或64位
np.int8|	8位的整形(-128~127)
np.int16|	-32768~32767
np.int32|	-2 ** 31 ~ 2 ** 31 - 1
np.int64|	-2 ** 63 ~ 2 ** 63 - 1
np.uint8|	8位的整形(0~255)
np.uint16|	-32768~32767
np.uint32	|0 ~ 2 ** 32 - 1
np.uint64|	0 ~ 2 ** 64 - 1
np.float16|	1位符号位，5位指数位，10位
np.float32	|1位符号位，8位指数位，23位
np.float64、np.float	|1位符号位，11位指数位，52位


#### 6.2.3 特定分布随机数创建矩阵函数

函数|	说明
--|--
np.random.rand(d0,d1,...,dn)   |	根据d0-dn(维度)创建随机数数组，[0,1)，均匀分布
np.random.randn(d0,d1,...,dn)	|根据d0-dn(维度)创建随机数数组，标准正态分布
np.random.randint(low[,high,shape])|	根据shape创建随机整数或整数数组，范围是[low,high)
np.random.normal(loc,scale,size)	|产生具有正态分布的数组，loc为均值，scale标准差，size为形状
np.random.permutation(a)	|根据数组a的第1轴产生一个新的乱序数组，不改变数组a
np.random.choice(a[,size,replace,p])|	从一维数组a中以概率p抽取元素，形成size形状新数组replace表示是否可能重用元素，默认为False
np.random.uniform(low,high,size)	|产生具有均匀分布的数组，low起始值，high结束值，size为形状
np.random.poisson(lam,size)|	产生具有泊松分布的数组，lam为随机事件发生率，size为形状

```
import numpy as np #为了方便使用numpy 采用np简写

#产生均匀分布的数据，2x3x4 三维数据
np.random.rand(2,3,4) 

# 产生20个均值为2、标准差为0.1满足正态分布的随机数序列
a = np.random.normal(2, 0.1, 20)

```



#### 6.2.4 矩阵/数组对象相关方法与属性 

##### reshape()  和resize()方法
- reshape()方法，在reshape方法里以元组、列表给出变化后的形状数据，并不影响原数组会新生成一个多维数组。
- resize() 方法，会修改数组本身的shape属性来改变数组的维度，原数组发生改变。

```
import numpy as np


#reshape方法，注意a被修改了
a = np.arange(12)
print(a)
a.reshape((3, 4))  
print(a)

#resize方法，注意a被修改了
a = np.arange(12)
a.resize([3, 4])
print(a)


#通过属性修改维度，注意a被修改了
a = np.arange(12)
print(a)
a.shape = (2, 6)


```

##### 数组的一些属性

```
# Array properties
import numpy as np
a = np.array([[11, 12, 13, 14, 15],
              [16, 17, 18, 19, 20],
              [21, 22, 23, 24, 25],
              [26, 27, 28 ,29, 30],
              [31, 32, 33, 34, 35]])

print(type(a)) # >>><class 'numpy.ndarray'>    正如你在上面的代码中看到的，NumPy数组实际上被称为：多维数组，ndarray。
print(a.dtype) # >>>int32    数据类型
print(a.size) # >>>25    数据总个数
print(a.shape) # >>>(5, 5)   数组的形状是它有多少行和列，上面的数组有5行和5列，所以它的形状是(5，5)。
print(a.ndim) # >>>2   数组的维数。

```



##### 矩阵数据支持索引，切片 类似列表

```
import numpy as np
a = np.arange(24)
print (a[3])

#二维数组索引
#对于2D数组：行的切片，列的切片。
a = np.arange(12)
a = a.reshape((3, 4))
print(a[0, 1:4])
print(a[1:4, 0])
print(a[::2,::2])

#三维数组索引
a = np.arange(24)
b = a.reshape((2, 3, 4))
print (b[1,1:2,1:])

```

###### 数组支持布尔筛选


```
import numpy as np
a = np.arange(0, 100, 10)
b = a[:5]
c = a[a >= 50]
print(b) 
print(c) 


# Where 函数筛选
a = np.arange(0, 100, 10)
b = np.where(a < 50) 
c = np.where(a >= 50)[0]
print(b) 
print(c) 
```



#### 6.2.5 NumPy的数学统计函数

函数|	说明
--|--
np.abs() np.fabs()|	计算数组各元素的绝对值
np.sqrt()|	计算数组各元素的平方根
np.square()	|计算数组各元素的平方
np.log(x),np.log10(x),np.log2(x)|	计算数组各元素的自然对数、10底对数和2底对数
np.ceil(x),np.floor(x)	|计算数组各元素的ceiling值或floor值
np.rint(x)|	计算数组各元素的四舍五入值
np.modf(x)	|将数据各元素的整数和小数部分以两个独立的数组形式返回
np.cos/cosh/sin/sinh/tan/tanh	|计算数据各元素的普通型和双典型的三角函数
np.exp(x)	|计算数组各元素的指数值
np.sum(a,axis=None)|	根据给定axis计算数组a相关元素之和，axis整数或元组
np.mean(a,axis=None)	|根据给定axis计算数组a相关元素的期望，axis整数或元组
np.average(a,axis=None,weights=None)|	根据给定axis计算数组a相关元素的加权平均值
np.std(a,axis=None)	|根据给定轴axis计算数组a相关元素的标准差
np.var(a,axis = None)	|根据给定轴axis计算数组a相关元素的方差
np.cov(a,axis = None)	|计算协方差
np.corrcoef	|计算相关系数，参阅
np.min(a) max(a)	|计算数组a中元素的最小值，最大值
np.argmin(a) argmax(a)	|计算数组a中元素的最小值，最大值的降一维后下标
np.unravel_index(index,shape)	|根据shape将一维下标index转换成多维下标
np.ptp(a)	|计算数组a中元素最大值和最小值的差
np.median(a)	|计算数组a中元素的中位数(中值)

iris.data测试数据下载：https://archive.ics.uci.edu/ml/machine-learning-databases/iris/
```
import numpy as np
import os
os.chdir("D://python_script//")
a, b = np.loadtxt("iris.data",  delimiter=',', usecols = [0, 1],unpack=True)
print(a)
print(b)

print(np.average(a))
print(np.mean(a))
print(np.max(a))
print(np.argmax(a))
print(np.min(a))
print(np.argmin(a), a[np.argmin(a)])
print(np.ptp(a))
print(np.std(a))
print(np.median(a))
print(np.sqrt(a), np.square(a))
print(np.rint(a))
print(np.ceil(a))


#print(np.column_stack([a, b]))
```
NumPy的算术运算，类似线性代数计算

```
import numpy as np
# Basic Operators
a = np.arange(25)
a = a.reshape((5, 5))

b = np.array([10, 62, 1, 14, 2, 56, 79, 2, 1, 45,
              4, 92, 5, 55, 63, 43, 35, 6, 53, 24,
              56, 3, 56, 44, 78])
b = b.reshape((5,5))

print(a + b)
print(a - b)
print(a * b)
print(a / b)
print(a ** 2)
print(a < b) 
print(a > b)
```


### 6.3 Pandas简介

基于NumPy 的一种工具，该工具是为了解决数据分析任务而创建的。Pandas 纳入了大量库和一些标准的数据模型，提供了高效地操作大型数据集所需的工具。最具有统计意味的工具包，某些方面优于R软件。数据结构有一维的Series，二维的DataFrame，三维的Panel。



- 官方网站：http://pandas.pydata.org/
- 中文文档：http://liao.cpython.org/partFive/

#### Pandas数据类型

- **Series**：一维数组，与Numpy中的一维array类似。二者与Python基本的数据结构List也很相近，其区别是：List中的元素可以是不同的数据类型，而Array和Series中则只允许存储相同的数据类型，这样可以更有效的使用内存，提高运算效率。
- **DataFrame**：二维的表格型数据结构。很多功能与R中的data.frame类似。可以将DataFrame理解为Series的容器。
- **Time-Series**：以时间为索引的Series。
- **Panel**：三维的数组，可以理解为DataFrame的容器。则可以视为excel的多表单sheet。




#### 6.3.1 安装
```
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple/  pandas

import pandas as pd
```

#### 6.3.2 创建Series

Series构造函数的name参数是给这列数据指定字段名。从结果可以看出t有两个名为'a'的label，值分别为2和7。
```
import pandas as pd
i = ["a", "c", "d", "a"]
v = [2, 4, 5, 7]
t = pd.Series(v, index=i, name = "col_name")
print(t)

#数据的访问，支持索引，切片，name取值
print("t[d]->", t["d"])
print('t[0 : 3]->', t[0 : 3])
print('t["a" : "c"]->', t["c" : "d"])  #注意label唯一才行，不然会报错

```

![image](https://note.youdao.com/yws/api/group/92990402/noteresource/9B633333F36D413FB1F1BEA6EC8709BC/version/2205?method=get-resource&shareToken=08DB9C16FB6D47AA8EBD104017423D2D&entryId=446741437)



#### Series对象的方法与属性

##### get方法
首先看看get方法，可以返回指定的key所对应的value值

```
import pandas as pd
idx =  "hello the cruel world".split()
val = [1, 21, 13, 104]
t = pd.Series(val, index = idx)

t.get("the")
t.get("The", "None")  #如果key不存在，返回default的值。

```

##### add、append方法
add和append方法都能改变series，只不过add类似于加法操作，而append则是连接。

```
import pandas as pd
idx =  "hello the cruel world".split()
val1 = [1, 21, 13, 104]

t = pd.Series(val1, index = idx)

val2 = [4, 4, 4, 4]
s = pd.Series(val2, index = idx)
t.add(s) 
t + 4

t.append(s)

```

#### 6.3.3 DataFrame创建

Pandas的Dataframe是二维的，每一列都是一个Series结构。


```
#手动创建
import pandas as pd
import numpy as np
df1 = pd.DataFrame({'A': 1.,
                    'B': pd.date_range('20130101', periods=4),
                    'C': pd.Series(1, index=list(range(4)), dtype='float32'),
                    'D': np.array([3] * 4, dtype='int32'),
                    'E': ["test", "train", "test", "train"],
                    'F': 'foo'},index=list("abcd"))
#numpy矩阵手动创建
df2 = pd.DataFrame(np.random.randn(10,3), columns = ["ca", "cb", "cc"], index =list("abcdefghij"))

```



#### DataFrame对象的属性与方法学习

属性|描述
--|--
columns |columns属性可以获得dataframe有那些列，即dataframe的index
shape |shape属性是描述dataframe的形状的
size  |dataframe的size属性返回的是dataframe的value的个数
values  |返回当前dataframe的数据和index、columns相对应。
dtypes |述当前dataframe的里的每列值的数据类型。
ndim |返回数据框维度
T  |dataframe的T属性，实际是转置的意思。

```
import pandas as pd
import numpy as np
df = pd.DataFrame(np.random.randn(10,3), columns = ["ca", "cb", "cc"], index =list("abcdefghij"))
df.columns
df.index
df.shape
df.size
df.values
df.dtypes
df.ndim
df.T

#方法
df.head()
df.tail(3)
df.to_numpy()
df.describe()  #方法显示数据的快速统计摘要
```


#### DataFrame数据访问与筛选
需要特别注意的是DataFrame和其他表格数据不一样的是，DataFrame是列访问机制。
```
import pandas as pd
import numpy as np
val = np.arange(10, 40).reshape(10, 3)
idx = ["ax", "bx", "cx"]
df = pd.DataFrame(val, columns = idx,index=list("abcdefghij"))

#[]号访问列
df["ax"]   #单列索引
df[["ax", "cx"]]  #多列手动选择

#利用列名访问列
df.ax

#DataFrame[start:end]则是通过切片选择的是行。
df["a" : "e"]
df[:3]    #多列切片


###利用方法进行选择

#按位置索引选择：iloc[]行列切片
df.iloc[1] #单独使用选择行
df.iloc[2 : 6, 0 : 2]   #行列选择
df.iloc[[0, 1, 3]]   #不同行选择
df.iloc[[0, 1, 3],[2,1]]   #不同行列选择


#按标签选择  loc[]行列切片， 行列的名字  dataFrame里可以通过loc[]的方式选择label标识的行数据。
df.loc["a"] #单独使用选择行
df.loc[["a","c"]]
df.loc[["a","c"],["ax","cx"]]
df.loc["b" : "e", "bx" : "cx"]
df.loc[: , "bx" : "cx"]


```

##### bool逻辑筛选数据

```
import pandas as pd
import numpy as np
val = np.arange(10, 60).reshape(10, 5)
col = ["ax", "bx", "cx", "dx", "ex"]
idx = list("abcdefghij")
df = pd.DataFrame(val, columns = col, index = idx)

#筛选，bx列大于30的数据
bs = df["bx"] > 30
df[bs]

#组合选择，bx列大于30的数据并且cx列大于40的数据

bs = (df["bx"] > 30) & (df["cx"] > 40)
df[bs]     #选择符合条件的行

#布尔选择的结果还是DataFrame，所以对于结果可以进行切片、label、loc等访问。

```


##### dataframe数据分类统计（重要）

apply 方法  和  groupby方法对数据进行分类统计
```
import pandas as pd
import numpy as np
val = np.arange(10, 60).reshape(10, 5)
col = ["ax", "bx", "cx", "dx", "ex"]
idx = list("abcdefghij")
df = pd.DataFrame(val, columns = col, index = idx)

df.apply(lambda col : col.sum(), axis = 0)  #求列的和
df.apply(lambda row : row.sum(), axis = 1)  #求行的和
df["plus"] = df.apply(lambda row : row.ax + row.cx, axis = 1) #两列相加，并增加新列，plus
df
```
groupby 方法，实现分类汇总
```
import pandas as pd
idx = [101,101,101,102,102,102,103,103,103]
name = ["apple","pearl","orange", "apple","pearl","orange","apple","pearl","orange"]
price = [1.0,2.0,3.0,4.00,5.0,6.0,7.0,8.0,9.0]
df0 = pd.DataFrame({ "fruit": name, "price" : price, "supplier" :idx})

dg =  df0.groupby("fruit")

#批量查看分组结果
for n, g in dg:
    print("group_name:", n, "\n",g,)

#分析结果总结
dg.describe()

#选择需要研究的列,属性或者方法获得统计结果
dg['price'].mean() 
dg['supplier'].value_counts()


#多个分组，大分组之后再分亚组
import pandas as pd
idx = [101,101,101,102,102,102,103,103,103,101,101,101,102,102,102,103,103,103]
name = ["apple","pearl","orange", "apple","pearl","orange","apple","pearl","orange","apple","pearl","orange", "apple","pearl","orange","apple","pearl","orange"]
price = np.arange(18)
df0 = pd.DataFrame({ "fruit": name, "price" : price, "supplier" :idx})

dg2 =  df0.groupby(["fruit", "supplier"])
for n, g in dg2:
    print("multiGroup on:", n, "\n",g)
dg2.describe()

```

==groupby方法总结==

>首先通过groupby得到DataFrameGroupBy对象, 然后选择需要研究的列,这样我们就得到了一个SeriesGroupby, 它代表每一个组都有一个Series
对SeriesGroupby进行操作, 比如.mean(), 相当于对每个组的Series求均值
#### 6.3.4 pandas 应用：读入写出数据处理数据

测试数据：https://archive.ics.uci.edu/ml/machine-learning-databases/iris/
```
import os
import pandas as pd
os.chdir("D://python_script//")
fn = "D://python_script//iris.data"
cols_name = ['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'class']
df = pd.read_csv(fn, names = cols_name)  #更多读取数据：https://pandas.pydata.org/pandas-docs/stable/reference/io.html

data_df=df.iloc[:,0:4]   #筛选一下数据

#利用apply计算：行列的和

data_df.apply(lambda col : col.sum(), axis = 0)
data_df.apply(lambda row : row.sum(), axis = 1)

#自定义计算
data_df["sepal_length x sepal_width"] =data_df.apply(lambda row : row.sepal_width * row.sepal_length, axis = 1)
data_df

data_df.to_csv("out.csv")  #计算结果写出

```



