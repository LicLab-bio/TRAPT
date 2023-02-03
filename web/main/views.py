
import main.models as md
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render


######## 网页路由 ########
def index(request):
    return render(request, 'index.html')

def analysis(request):
    return render(request, 'analysis.html')

def analysis_trs_activity(request):
    return render(request, 'analysis_trs_activity.html')

def download(request):
    return render(request, 'download.html')

def help(request):
    return render(request, 'help.html')

##### 进程路由 #####
def isFinal(request):
    uuid = request.POST.get('uuid')
    result = md._isFinal(uuid=uuid)
    return JsonResponse(result, safe=False)

def run(request):
    gene_set = request.POST.get('gene_set')
    vgae_iterations = request.POST.get('vgae_iterations')
    range_threshold = request.POST.get('range_threshold')
    uuid = md._run(gene_set=gene_set,vgae_iterations=vgae_iterations,range_threshold=range_threshold)
    return JsonResponse(uuid, safe=False)

def tryAgain(request):
    uuid = request.POST.get('uuid')
    result = md._tryAgain(uuid=uuid)
    return JsonResponse(result, safe=False)


######## 数据路由 ########

# TRs Activity Table 数据
def getTRsActivityTable(request):
    uuid = request.POST.get('uuid')
    result = md._getTRsActivityTable(uuid=uuid)
    return JsonResponse(result, safe=False)

# 运行日志
def getLog(request):
    uuid = request.POST.get('uuid')
    result = md._getLog(uuid=uuid)
    return HttpResponse(result)

# jbrowse
def drawLine(request):
    uuid = request.POST.get('uuid')
    tr = request.POST.get('tr')
    result = md._drawLine(uuid=uuid,tr=tr)
    return JsonResponse(result)
