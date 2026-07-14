Set-Location -LiteralPath $PSScriptRoot
if (!(Test-Path -LiteralPath "logs")) {
    New-Item -ItemType Directory -Path "logs" | Out-Null
}
matlab -batch "cd('$PSScriptRoot'); run_HSM_runtime_ablation" -logfile "logs\HSM_analysis_launch.log"
