from django.urls import path

from main import views

urlpatterns = [
    ##### 网页路由 #####
    path('TRAPT', views.index, name='index'),
    path('TRAPT/analysis', views.analysis, name='analysis'),
    path('TRAPT/browse', views.browse, name='browse'),
    path('TRAPT/analysis_trs_activity', views.analysis_trs_activity, name='analysis_trs_activity'),
    path('TRAPT/download', views.download, name='download'),
    path('TRAPT/help', views.help, name='help'),
    ##### TRAPT/进程路由 #####
    path('TRAPT/isFinal', views.isFinal, name='isFinal'),
    path('TRAPT/run', views.run, name='run'),
    path('TRAPT/tryAgain', views.tryAgain, name='tryAgain'),
    ##### TRAPT/数据路由 #####
    path('TRAPT/getLog', views.getLog, name='getLog'),
    path('TRAPT/drawLine', views.drawLine, name='drawLine'),
    path('TRAPT/getTRsInfo', views.getTRsInfo, name='getTRsInfo'),
    path('TRAPT/getTRsActivityTable', views.getTRsActivityTable, name='getTRsActivityTable'),
]
