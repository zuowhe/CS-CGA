Set-Location -LiteralPath "E:\GPT projects\BNSL_CSCGA"
if (!(Test-Path -LiteralPath "logs")) {
    New-Item -ItemType Directory -Path "logs" | Out-Null
}
matlab -batch "cd('E:\GPT projects\BNSL_CSCGA'); Main_HSM_Analysis" -logfile "logs\HSM_analysis_launch.log"
