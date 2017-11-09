from django import forms

from .models import single_loci_blast
from .models import multiple_loci_blast

class SLB_form(forms.ModelForm):

    class Meta:
        model = single_loci_blast
        fields = ('sample_name','loci_name','seq',)

class MLB_form(forms.ModelForm):

    class Meta:
        model = multiple_loci_blast
        fields = ('seq',)