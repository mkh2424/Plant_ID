from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name='home'),
    url(r'^index/$', views.index),
    url(r'^SLB/$',views.single_loci_blast),
    url(r'^MLB/$',views.multiple_loci_blast),
    url(r'^about/$',views.about)
]