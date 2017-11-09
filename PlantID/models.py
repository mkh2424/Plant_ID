from django.db import models
from django.utils import timezone

Barcode_loci_choice = (("accD", "accD"), ("matK", "matK"), ("ndhJ", "ndhJ"), ("rpoB","rpoB"),
                ("rpoC1", "rpoC1"), ("ycf5","ycf5"), ("rbcL-a", "rbcL-a"), ("trnH-psbA", "trnH-psbA"),
                ("ITS1", "ITS1"), ("ITS2", "ITS2"))

class single_loci_blast(models.Model):
    sample_name = models.CharField(max_length=100)
    loci_name = models.CharField(max_length= 15, choices=Barcode_loci_choice,default="accD")
    seq = models.TextField()

class multiple_loci_blast(models.Model):
    seq = models.TextField()

# Create your models here.
