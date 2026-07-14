@echo off
cd /d "%~dp0"
if not exist logs mkdir logs
matlab -batch "cd('E:\GPT projects\BNSL_CSCGA'); Main_HSM_Analysis" > "logs\HSM_analysis_launch.log" 2>&1
