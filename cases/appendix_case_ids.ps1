# Shared case identity parsing. Preserve published APP identifiers when headings move.
function Convert-AppendixVertex([string]$Text) {
    $vector = @(0,0,0,0,0)
    $names = @('aaa','bbb','ccc','ddd','eee')
    $hits = [regex]::Matches($Text, '\\(aaa|bbb|ccc|ddd|eee)(?:\^\{?([0-9]+)\}?)?')
    if ($hits.Count -eq 0) { throw "Unparsed vertex '$Text'" }
    foreach ($hit in $hits) {
        $index = [array]::IndexOf($names, $hit.Groups[1].Value)
        $vector[$index] += $(if ($hit.Groups[2].Success) { [int]$hit.Groups[2].Value } else { 1 })
    }
    return '[' + ($vector -join ',') + ']'
}

function Get-AppendixVertexKey([string]$Vectors) {
    $parts = @([regex]::Matches($Vectors, '\[[0-9]+(?:,[0-9]+){4}\]') | ForEach-Object { $_.Value })
    if ($parts.Count -ne 3) { throw "Expected three vertex vectors: $Vectors" }
    return ($parts | Sort-Object) -join ';'
}

function Read-AppendixCaseIds([string]$Path) {
    $ids = @{}
    $seenIds = @{}
    foreach ($line in (Get-Content -LiteralPath $Path)) {
        $match = [regex]::Match($line, '^\s*\["(APP-[^"]+)",\s*(\[\[[0-9,]+\],\[[0-9,]+\],\[[0-9,]+\]\])')
        if (-not $match.Success) { continue }
        $id = $match.Groups[1].Value
        $key = Get-AppendixVertexKey $match.Groups[2].Value
        if ($ids.ContainsKey($key) -or $seenIds.ContainsKey($id)) {
            throw "Duplicate identity in reference matrix: $id"
        }
        $ids[$key] = $id
        $seenIds[$id] = $true
    }
    if ($ids.Count -eq 0) { throw "No APP identities in $Path" }
    return $ids
}

function Get-AppendixHeading([string]$Line, [hashtable]$Ids) {
    if ($Line -notmatch '^\s*\\subsubsection\*\{Case\.?\s') { return $null }
    $match = [regex]::Match($Line, 'Case\.?\s+\$\\\{(.*?)\\\}')
    if (-not $match.Success) { throw 'No heading boundary' }
    $body = $match.Groups[1].Value
    $primary = ($body.Replace('\,','') -split '\(\\text\{or\}')[0]
    $parts = @($primary -split ',' | Where-Object { $_.Trim().Length -gt 0 })
    if ($parts.Count -eq 2) { return [pscustomobject]@{ VertexCount=2; Id=$null; Body=$body } }
    if ($parts.Count -ne 3) { throw "Primary vertex count=$($parts.Count)" }
    $vectors = @($parts | ForEach-Object { Convert-AppendixVertex $_ })
    $key = Get-AppendixVertexKey ('[' + ($vectors -join ',') + ']')
    if (-not $Ids.ContainsKey($key)) { throw "Unknown case; assign a stable ID in the reference matrix: $body" }
    $id = $Ids[$key]
    if ($id -match '-OR') { throw "An OR branch became a primary heading: $id" }
    return [pscustomobject]@{ VertexCount=3; Id=$id; Body=$body }
}
