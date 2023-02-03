from django.urls import path
from main import views

urlpatterns = [
    ##### 网页路由 #####
    path('', views.index, name='index'),
    path('analysis', views.analysis, name='analysis'),
    path('analysis_trs_activity', views.analysis_trs_activity, name='analysis_trs_activity'),
    path('download', views.download, name='download'),
    path('help', views.help, name='help'),
    ##### 进程路由 #####
    path('isFinal', views.isFinal, name='isFinal'),
    path('run', views.run, name='run'),
    path('tryAgain', views.tryAgain, name='tryAgain'),
    ##### 数据路由 #####
    path('getLog', views.getLog, name='getLog'),
    path('drawLine', views.drawLine, name='drawLine'),
    path('getTRsActivityTable', views.getTRsActivityTable, name='getTRsActivityTable'),
]
