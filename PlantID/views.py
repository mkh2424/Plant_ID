from django.shortcuts import render
from django.http import HttpResponse
from .forms import SLB_form
from .forms import MLB_form
from .multipe_loci_blast import *
from .single_loci_blast import *
from os import remove

from .models import single_loci_blast
from .models import multiple_loci_blast


# Create your views here.

def index(request):
    return render(request, 'PlantID/index.html', {})

def single_loci_blast(request):
    if request.method == "POST":
        form = SLB_form(request.POST)
        if form.is_valid():
            single_loci_blast_wrapper(form.cleaned_data['sample_name'], form.cleaned_data['loci_name'], form.cleaned_data['seq'])
            blastn = ""
            blastx = ""
            with open("Output_"+form.cleaned_data['sample_name']+"_blastn.xml", 'r') as data_handle:
                blastn = data_handle.read().replace('\n','<br />')
            with open("Output_"+form.cleaned_data['sample_name']+"_blastx.xml", 'r') as data_handle:
                blastx = data_handle.read().replace('\n','<br />')
            return HttpResponse(blastn + "<br />" + blastx)
    else:
        form = SLB_form()
    return render(request, 'PlantID/single_loci_blast.html',{'form' : form})

def multiple_loci_blast(request):
    if request.method == "POST":
        form = MLB_form(request.POST)
        if form.is_valid():
            multiple_loci_blast_wrapper(form.cleaned_data['seq'])
            with open("output_test.txt", 'r') as data_handle:
                #remove("input_test.fasta")
                #remove("output_test.txt")
                #remove("parsed_test.fasta")
                return HttpResponse(data_handle.read().replace('\n','<br />'))
    else:
        form = MLB_form()
    return render(request, 'PlantID/multiple_loci_blast.html', {'form': form})

def about(request):
    return render(request,'PlantID/about.html')