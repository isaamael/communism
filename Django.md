# Django 的使用全部基于python、pycharm 和Django环境



### 0.Django的安装通过pycharm pro版完成



### 1.建立Django project

#通过pycharm pro 直接创建DjangoProject1

#创建一个app，命名为DjangoWeb

```powershell
python manage.py startapp DjangoWeb 
```

#修改DjangoWeb这个app中的views.py，设置一个测试函数

```python
from django.http import HttpResponse
from django.shortcuts import render

# Create your views here.
def index(request):
    return HttpResponse('I get it')
```

#在根目录的urls.py 中添加views.index这个方法

```python
from DjangoWeb import views

urlpatterns = [
    path("admin/", admin.site.urls),
    path("index/", views.index),
]
```

#修改根目录下settings.py，把DjangoWeb这个app挂载到INSTALLED_APPS下

```python
INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "DjangoWeb",
]
```

#powershell 运行 启动Django的本地维护

```powershell
python manage.py runserver
```

#修改根目录下settiing.py中的语言源码，界面显示成中文

```python
LANGUAGE_CODE = "zh-hans"
```

#在浏览器的localhost:8000/index下可以显示views.index的内容，完成基本app配置



#添加超级管理员账户，方便后台管理

```powershell
python manage.py createsuperuser
```

#输入用户名和密码（都是isaac）

#在地址localhost:8000/admin 下登录

#进入网站后台

#以上完成DjangoProject1的初步准备，python manage.py help 可以查询命令信息



### 2.添加article

#article的添加通过DjangoWeb下的model.py完成，修改该文件

```python
from django.db import models

# Create your models here.
#添加类
class article(models.Model):   #创建article 类，继承父类，包括CharField和TextField
    title=models.CharField(max_length=100) #为article 添加属性title
    content=models.TextField()  #为article 添加属性content
```

#在model.py中添加类以后，需要将这个模型迁移到数据库中，通过两步完成

```powershell
python manage.py makemigration   #制造迁移

python manage.py migrate   #迁移
```

#此时，article模型已经被迁移到数据库中，数据库在根目录下的db.sqlite3

#要在后台地址看到article表，仍需要在DjangoWeb这个app的admin.py中注册article 模型

```python
from django.contrib import admin
from DjangoWeb.models import article  #加载DjangoWeb路径下model模块的article类
# Register your models here.
admin.site.register(article)   #注册这个类，使之可以通过admin后台管理
######看着这个 类 熟悉的界面，原来期刊们的投稿主页都是用Django做的= = ##############
```



#此时article 类虽然已经部署，需要通过一个处理方法解决所有article的实例，通过每一个实例的唯一标识 id 来解决

```python
#模型的objects是获取或操作模型的对象
#article.objects.get(条件)
#article.objects.all()
#article.objects.filter(条件)
```

#在views.py中定义新的处理article的函数

```python
from django.http import HttpResponse,Http404
from django.shortcuts import render
from DjangoWeb.models import article
# Create your views here.
def index(request):
    return HttpResponse('I get it')
#定义处理函数
def article_detail(request, article_id):   #添加参数，article id
    global article1
    try:
        article1=article.objects.get(id=article_id)  #通过object方法，获得article 类 中，id与article id一致的实例，赋予article1
    except article1.DoesNotExist:   #如果无法取出实例，就调用404
        raise Http404("not exist")
    #成功取出实例，返回表达式的结果，article1为实例，具有module中定义的title和content属性
    return HttpResponse("<h2>文章标题: %s </h2><br>文章内容： %s" % (article1.title, article1.content))  
```

 #同时在路由器urls.py中定义变量访问article的路径

```python
from DjangoWeb import views

urlpatterns = [
    path("admin/", admin.site.urls),
    path("index/", views.index),
    #访问article实例的路径是变量
    #同时通过将客户端访问时的article id 传为 view 中定义的 article_detail函数中的 article_id参数，实现调用相应实例
    path("article/<int:article_id>",views.article_detail,name="article_detail")
]
```

#整体逻辑就是将唯一识别的id作为变量符号，在客户端访问路径和处理函数中保持一致，从而实现动态访问所有实例

#有些期刊的网站也是这个方法，也有一些是通过doi号码，应该是规定的适合，将id替换为了doi，毕竟二者同样都是一个article实例 唯一识别的信息



#之后进行前端后端的分离，及路由分发的分离

#这部分内容需要理清每一个py文件所接收和指向的位置

#Django直接实现的前后端分离，是基于将 ‘处理方法’ 转嫁到 ‘模板’ 文件夹中的前端渲染实现的

#在 ‘app’ 下新建 ‘templates’ 文件夹，在下面建每一个页面对应的 html 文件，此时该文件即为一个前端渲染文件

#在 html 文件中接收 ‘view’ 中传递的内容并进行 ‘加工显示’ ，这部分内容通过 http 语法实现

```python
#app中 view 的内容

from django.http import HttpResponse, Http404
from django.shortcuts import render, get_object_or_404
from DjangoWeb.models import article

# Create your views here.
#定义index显示的内容
def index(request):
    return HttpResponse('I get it')

#定义处理方法，获得文章内容
def article_detail(request, article_id):
    article1=get_object_or_404(article, pk=article_id)    #获得每个文章的内容
    context = {}
    context["article_obj"] = article1    #将文章内容写入字典
    return render (request, "article.detail.html", context)   #render函数把字典返还到一个负责前端渲染的 html 中

#定义处理方法，获得文章列表
def article_list(request):
    articles=article.objects.all()    # objects.all 获得文章全部信息
    context={}
    context["articles"]=articles   #写入到字典中
    return render(request,"article.list.html",context)   #同样通过render函数返还到前端的 html 中

#通过以上方法，收到request后，把相应的内容传递到相应的 html 中
```

```http
#显示每个文章的标题和内容的html
#article.detail.html的内容
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Title</title>
</head>
<body>
    #写出相应的渲染方法，如何显示所传递的内容
    <h2>{{ article_obj.title }}</h2> 
    <hr>
    <p>{{ article_obj.content }}</p>
</body>
</html>
```

```http
#显示每个文章的标题和内容的html
#article.list.html的内容
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>article</title>
</head>
<body>
{% for article in articles %}
    #这里的渲染既显示了需要的内容，同时通过http中的<a href>语句，对路由进行了分发
    <a href={% url "article_detail" article.pk  %}>{{ article.title }}</a>
{% endfor %}
</body>
</html>
```

#路由的分发通过类似子文件夹的方法实现，即 在 全局的 urls.py 中 挂载一个 子app中的 urls.py，保证流入和指出 畅通即可

#实现结果如下

```python
#全局的urls
from django.contrib import admin
from django.urls import path, include
from DjangoWeb import views

urlpatterns = [
    path("admin/", admin.site.urls),
    path("index/", views.index),
    path("article/",include("DjangoWeb.urls")) #这里进行了挂载子路径，挂载通过 include 方法
]
```

```python
#app中的子路由
from django.urls import path
from DjangoWeb import views

urlpatterns = [
    #子路由的路径与全局路由的路径是衔接的
    #下面这个path实际的绝对路径是 localhost:8000/article/<int:article_id>
    path("<int:article_id>",views.article_detail,name="article_detail"), 
    path("", views.article_list, name="article_list")
]
```



#以上内容

#通过将view.py中的处理和渲染两部分内容分开，将渲染传递到单独的渲染端，实现前后端的分离

#通过渲染时的 http 语句，实现显示内容和路由地址的耦联

#通过将总路由与app内的子路由连接在一起，实现路由的分发









