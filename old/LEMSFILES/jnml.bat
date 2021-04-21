@echo off

REM
REM  This file can be used to run the jNeuroML jar on Windows. Type:
REM   
REM     jnml.bat
REM   
REM  inside this directory, or add the path to this directory to your PATH 
REM  environment variable and create the variable JNML_HOME, pointing to this 
REM  directory, so that you can run: 
REM  
REM      jnml.bat
REM  
REM  from any directory
REM


REM  Set the current version of jNeuroML
set VERSION=0.10.1


REM Create the Java classpath
set CLASSPATH=jNeuroML-%VERSION%-jar-with-dependencies.jar;%JNML_HOME%\jNeuroML-%VERSION%-jar-with-dependencies.jar


REM Run this file with the java executable, passing on arguments given to jnml
java -Xmx400M -cp %CLASSPATH% org.neuroml.JNeuroML %1 %2 %3 %4 %5 %6
