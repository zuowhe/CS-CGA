@echo off
cd /d "%~dp0"
if not exist logs mkdir logs
matlab -batch "cd('%~dp0'); run_HSM_runtime_ablation" > "logs\HSM_analysis_launch.log" 2>&1
