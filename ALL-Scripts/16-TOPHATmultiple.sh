#!/bin/sh
#PBS -N toPHATPIPELINE
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=300:00:00,cput=300:00:00,nodes=1
#PBS -q long
#PBS -d /home/amodupe/CARL/TORQUEoutput

time perl /home/amodupe/CARL/Scripts/15-tophatmultiple.pl
