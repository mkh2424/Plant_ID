from django.contrib import admin
from .models import single_loci_blast
from .models import multiple_loci_blast

admin.site.register(single_loci_blast)
admin.site.register(multiple_loci_blast)

# Register your models here.
