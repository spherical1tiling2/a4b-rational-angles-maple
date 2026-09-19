param(
    [string]$Maple = 'C:/Program Files/Maple 2024/bin.X86_64_WINDOWS/cmaple.exe',
    [string]$OutputDirectory = (Join-Path $PSScriptRoot 'results-new'),
    [string]$ManuscriptPath,
    [ValidateRange(1,16)][int]$Workers = 4,
    [ValidateRange(1,257)][int]$BatchSize = 16
)
$ErrorActionPreference = 'Stop'
$source = Join-Path $PSScriptRoot 'cases'
$output = [IO.Path]::GetFullPath($OutputDirectory)
$work = $output + '.work'
if ((Test-Path -LiteralPath $output) -or (Test-Path -LiteralPath $work)) {
    throw 'Choose a new output directory. Existing computations are never reused.'
}
if (-not (Test-Path -LiteralPath $Maple -PathType Leaf)) { throw 'Set -Maple to the cmaple executable.' }
New-Item -ItemType Directory -Path $output,$work | Out-Null
$started = [DateTime]::UtcNow
$capacity = [long](Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory
$drive = Get-PSDrive -Name ([IO.Path]::GetPathRoot($output).Substring(0,1))
$limits = [ordered]@{
    memory_capacity_bytes=$capacity
    soft_memory_bytes=[long](0.75*$capacity)
    hard_memory_bytes=[long](1.5*$capacity)
    disk_reserve_bytes=[long](0.05*($drive.Used+$drive.Free))
    limit_source='host_capacity'
}
$limits | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $output 'resource_limits.json') -Encoding utf8
$inputHashes = @(Get-ChildItem -LiteralPath $source -File | Sort-Object Name | ForEach-Object {
    [pscustomobject]@{path=('cases/'+$_.Name);sha256=(Get-FileHash -LiteralPath $_.FullName).Hash.ToLowerInvariant()}
})
$inputHashes | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $output 'input_hashes.json') -Encoding utf8
if ($ManuscriptPath) {
    (Get-FileHash -LiteralPath $ManuscriptPath).Hash.ToLowerInvariant() | Set-Content -LiteralPath (Join-Path $output 'manuscript.sha256') -Encoding ascii
    & (Join-Path $source 'test_case_extraction.ps1') -TexPath $ManuscriptPath | Set-Content -LiteralPath (Join-Path $output 'extraction_tests.log') -Encoding utf8
    & (Join-Path $source 'extract_appendix_case_matrix.ps1') -TexPath $ManuscriptPath -MatrixPath (Join-Path $work 'matrix.txt') -IssuePath (Join-Path $work 'matrix_issues.txt') | Out-Null
    & (Join-Path $source 'extract_appendix_reference_polynomials.ps1') -TexPath $ManuscriptPath -OutputPath (Join-Path $work 'polynomials.mpl') -IssuePath (Join-Path $work 'polynomial_issues.txt') | Out-Null
    & (Join-Path $source 'extract_appendix_solutions.ps1') -TexPath $ManuscriptPath -OutputPath (Join-Path $work 'candidates.txt') | Out-Null
    foreach ($pair in @(@('matrix.txt','appendix_all_3vertex_cases.txt'),@('polynomials.mpl','appendix_reference_polys_by_heading.mpl'),@('candidates.txt','appendix_3vertex_solution_candidates.txt'))) {
        $generated = (Get-Content -Raw -LiteralPath (Join-Path $work $pair[0])) -replace '\r\n',"`n"
        $included = (Get-Content -Raw -LiteralPath (Join-Path $source $pair[1])) -replace '\r\n',"`n"
        if ($generated.Trim() -ne $included.Trim()) { throw "Manuscript input mismatch: $($pair[1])" }
    }
}
$pending = [System.Collections.Generic.Queue[object]]::new()
$active = [System.Collections.Generic.List[object]]::new()
$completed = [System.Collections.Generic.List[object]]::new()
$peak = 0L
$softCrossed = $false
for ($first=1; $first -le 257; $first+=$BatchSize) {
    $last=[math]::Min(257,$first+$BatchSize-1)
    $name='batch-{0:D3}-{1:D3}' -f $first,$last
    $directory=Join-Path $work $name
    New-Item -ItemType Directory -Path $directory | Out-Null
    Get-ChildItem -LiteralPath $source -File | Copy-Item -Destination $directory
    Set-Content -LiteralPath (Join-Path $directory 'solution_range.mpl') -Value @("first_case := ${first}:","last_case := ${last}:") -Encoding ascii
    $pending.Enqueue([pscustomobject]@{first=$first;last=$last;name=$name;directory=$directory;process=$null;started=$null})
}
try {
    while ($pending.Count -or $active.Count) {
        while ($pending.Count -and $active.Count -lt $Workers -and -not $softCrossed) {
            $job=$pending.Dequeue(); $job.started=[DateTime]::UtcNow
            $job.process=Start-Process -FilePath $Maple -ArgumentList @('-q','run_solution_search.mpl') -WorkingDirectory $job.directory -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $job.directory 'stdout.log') -RedirectStandardError (Join-Path $job.directory 'stderr.log')
            $active.Add($job)
        }
        $processes=@(Get-CimInstance Win32_Process)
        $descendants=[System.Collections.Generic.HashSet[int]]::new()
        foreach ($job in $active) { [void]$descendants.Add($job.process.Id) }
        do {
            $added=$false
            foreach ($proc in $processes) {
                if ($descendants.Contains([int]$proc.ParentProcessId) -and $descendants.Add([int]$proc.ProcessId)) { $added=$true }
            }
        } while ($added)
        $rss=[long](($processes | Where-Object {$descendants.Contains([int]$_.ProcessId)} | Measure-Object WorkingSetSize -Sum).Sum)
        $peak=[math]::Max($peak,$rss)
        $softCrossed=$rss -ge $limits.soft_memory_bytes
        $free=(Get-PSDrive -Name $drive.Name).Free
        if ($rss -ge $limits.hard_memory_bytes -or $free -lt $limits.disk_reserve_bytes) { throw 'Run resource limit exceeded' }
        foreach ($job in @($active.ToArray())) {
            $job.process.Refresh()
            if (-not $job.process.HasExited) { continue }
            $stdout=Get-Content -Raw -LiteralPath (Join-Path $job.directory 'stdout.log')
            $stderr=Get-Content -Raw -LiteralPath (Join-Path $job.directory 'stderr.log')
            $m=[regex]::Match($stdout,'SUMMARY processed=(\d+) failed=(\d+) appendix=(\d+) found=(\d+) unique_found=(\d+)')
            if ($job.process.ExitCode -ne 0 -or -not $m.Success -or [int]$m.Groups[1].Value -ne ($job.last-$job.first+1) -or $m.Groups[2].Value -ne '0' -or ($stdout+$stderr) -match '(?im)^Error|syntax error|FAIL:|CASE_ERROR') {
                throw "Batch $($job.name) failed; inspect its logs"
            }
            $target=Join-Path $output ('search_batches/'+$job.name)
            New-Item -ItemType Directory -Path $target -Force | Out-Null
            foreach ($file in @('stdout.log','stderr.log','solution_candidates.csv','degree_bound_exclusions.csv','solution_factors.csv','solution_progress.log','solution_checkpoint.txt')) {
                if (Test-Path -LiteralPath (Join-Path $job.directory $file)) { Copy-Item -LiteralPath (Join-Path $job.directory $file) -Destination $target }
            }
            $completed.Add([pscustomobject]@{first=$job.first;last=$job.last;started_utc=$job.started.ToString('o');finished_utc=[DateTime]::UtcNow.ToString('o');seconds=([DateTime]::UtcNow-$job.started).TotalSeconds;processed=[int]$m.Groups[1].Value;found=[int]$m.Groups[4].Value})
            [void]$active.Remove($job)
            Write-Output "Completed $($job.first)-$($job.last); batches=$($completed.Count)"
        }
        $state=[ordered]@{status='running';started_utc=$started.ToString('o');updated_utc=[DateTime]::UtcNow.ToString('o');elapsed_seconds=([DateTime]::UtcNow-$started).TotalSeconds;peak_rss_bytes=$peak;measurement_scope='search_process_trees';workers=$Workers;completed=@($completed.ToArray())}
        $state | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $output 'execution.json') -Encoding utf8
        if ($active.Count) { Start-Sleep -Seconds 2 }
    }
    $verify=Join-Path $work 'verification'
    New-Item -ItemType Directory -Path $verify | Out-Null
    Get-ChildItem -LiteralPath $source -File | Copy-Item -Destination $verify
    $found=@(Get-ChildItem -Path (Join-Path $output 'search_batches/*/solution_candidates.csv') | ForEach-Object {Import-Csv -LiteralPath $_.FullName})
    $excluded=@(Get-ChildItem -Path (Join-Path $output 'search_batches/*/degree_bound_exclusions.csv') | ForEach-Object {Import-Csv -LiteralPath $_.FullName})
    $found | Sort-Object case,f | Export-Csv -LiteralPath (Join-Path $output 'search_candidates.csv') -NoTypeInformation -Encoding utf8
    $excluded | Sort-Object case,f | Export-Csv -LiteralPath (Join-Path $output 'degree_bound_exclusions.csv') -NoTypeInformation -Encoding utf8
    $mapleData=@('# Candidates and degree-bound exclusions from this execution.')
    foreach ($item in @(@('search_candidates',$found),@('degree_excluded_candidates',$excluded))) {
        $rows=@($item[1] | ForEach-Object {'["'+$_.case+'",'+$_.f+','+$_.angles_pi+']'})
        $mapleData+=($item[0]+" := [`n"+($rows -join ",`n")+"`n]:")
    }
    [IO.File]::WriteAllText((Join-Path $verify 'search_candidates.mpl'),($mapleData -join "`n")+"`n",[Text.UTF8Encoding]::new($false))
    Copy-Item -LiteralPath (Join-Path $verify 'search_candidates.mpl') -Destination $output
    foreach ($entry in @(@('validate_case_manifest.mpl','PASS: all affine equations'),@('run_polynomial_comparison.mpl','SUMMARY squarefree_same=230 shared=0 different=0 phase_deficient=0'),@('run_parameter_families.mpl','PASS: parameter_family_records=3'),@('run_verified_comparison.mpl','PASS: all finite manuscript records agree'),@('run_two_vertex_case.mpl','PASS: resultant_branches=15; exact_candidates=4; retained_manuscript_candidates=1'))) {
        $log=Join-Path $output ($entry[0].Replace('.mpl','.log'))
        Push-Location $verify
        try { $text=(& $Maple -q $entry[0] 2>&1 | Out-String); $code=$LASTEXITCODE } finally { Pop-Location }
        [IO.File]::WriteAllText($log,$text,[Text.UTF8Encoding]::new($false))
        if ($code -ne 0 -or -not $text.Contains($entry[1]) -or $text -match '(?im)^Error|syntax error|FAIL:') { throw "Verification failed: $($entry[0])" }
        Write-Output $text.Trim()
    }
    foreach ($file in @('appendix_polynomial_squarefree_comparison.csv','parameter_families.csv','family_lift_candidates.mpl','finite_candidate_checks.csv','manuscript_finite_candidates.csv','two_vertex_resultants.mpl','two_vertex_candidates.csv')) {
        Copy-Item -LiteralPath (Join-Path $verify $file) -Destination $output
    }
    $coverage=@(Get-ChildItem -Path (Join-Path $output 'search_batches/*/solution_progress.log') | ForEach-Object {
        $batch=$_.Directory.Name
        foreach ($line in (Get-Content -LiteralPath $_.FullName)) {
            if ($line -match '^CASE_BEGIN (APP-\d+(?:-OR\d+)?) raw_index=(\d+)') {
                [pscustomobject]@{case=$Matches[1];input_index=[int]$Matches[2];batch=$batch;status='complete'}
            }
        }
    })
    if ($coverage.Count -ne 257 -or @($coverage.case | Sort-Object -Unique).Count -ne 257) { throw 'Incorrect case coverage' }
    $coverage | Sort-Object input_index | Export-Csv -LiteralPath (Join-Path $output 'case_coverage.csv') -NoTypeInformation -Encoding utf8
    $state.status='complete'; $state.updated_utc=[DateTime]::UtcNow.ToString('o'); $state.elapsed_seconds=([DateTime]::UtcNow-$started).TotalSeconds
    $state.search_candidates=$found.Count; $state.degree_exclusions=$excluded.Count
    $state | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $output 'execution.json') -Encoding utf8
    Write-Output 'COMPLETE: 257 three-vertex cases and one two-vertex case; all comparisons passed.'
} catch {
    if ($descendants) { foreach ($processIdToStop in $descendants) { Stop-Process -Id $processIdToStop -Force -ErrorAction SilentlyContinue } }
    throw
}
