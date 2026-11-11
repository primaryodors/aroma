<?php
// Share the receptor page's empirical ordering with the docking entry point.
require_once(__DIR__ . '/protutils.php');
require_once(__DIR__ . '/odorutils.php');

$receptor = $argv[1] ?? '';
$ordered = array_keys(all_empirical_pairs_for_receptor($receptor));
$top = $prots[$receptor]['best_agonist'] ?? null;
if ($top && isset($odors[$top])) array_unshift($ordered, $top);

// The page appends prediction-only ligands after the empirical rows.
$prediction_file = __DIR__ . '/../predict/dock_results.json';
$predictions = is_file($prediction_file) ? json_decode(file_get_contents($prediction_file), true) : [];
$scores = [];
foreach (($predictions[$receptor] ?? []) as $name => $dock)
{
    if (!isset($dock['DockScore'])) continue;
    $odor = find_odorant($name);
    if ($odor) $scores[$odor['oid']] = floatval($dock['DockScore']);
}
arsort($scores);
$ordered = array_merge($ordered, array_keys($scores));
echo json_encode(array_values(array_unique($ordered))) . "\n";
