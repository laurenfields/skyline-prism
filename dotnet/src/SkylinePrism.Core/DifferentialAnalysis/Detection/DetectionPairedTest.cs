using System;
using System.Collections.Generic;
using System.Linq;

namespace SkylinePrism.Core.DifferentialAnalysis.Detection;

/// <summary>
/// One peptide's paired detection result.
/// </summary>
/// <param name="PeptideId">The peptide.</param>
/// <param name="Pairs">Matched subjects contributing.</param>
/// <param name="DetA">Pairs where it was detected in arm A.</param>
/// <param name="DetB">Pairs where it was detected in arm B.</param>
/// <param name="OnlyA">Discordant pairs detected in A but not B.</param>
/// <param name="OnlyB">Discordant pairs detected in B but not A.</param>
/// <param name="RateA">Detection rate in arm A.</param>
/// <param name="RateB">Detection rate in arm B.</param>
/// <param name="P">Two-sided exact McNemar p.</param>
/// <param name="Q">Adjusted p.</param>
public sealed record DetectionPairedRow(
    string PeptideId, int Pairs, int DetA, int DetB, int OnlyA, int OnlyB,
    double RateA, double RateB, double P, double Q);

/// <summary>
/// Paired detection: whether a peptide is detected at different rates between two conditions
/// measured on the same subjects.
/// </summary>
/// <remarks>
/// The paired counterpart of <see cref="DetectionTest"/>. Using the unpaired Fisher test on a paired
/// design throws away the pairing, which is the whole of the design's power - and, where subjects
/// differ a lot from each other, most of the signal with it. Only the discordant pairs contribute,
/// which is also why a paired detection result can rest on far fewer observations than its subject
/// count suggests; <see cref="DetectionPairedRow.OnlyA"/> and <see cref="DetectionPairedRow.OnlyB"/>
/// are reported so that is visible rather than implied.
/// </remarks>
// Public, like DetectionTest and DetectionGlm beside it: the GUI calls all three directly, and a
// detection test the pane cannot reach is not a detection test.
public static class DetectionPairedTest
{
    /// <summary>
    /// Run McNemar's test per peptide over the matched pairs, then adjust.
    /// </summary>
    /// <param name="detectionMatrix">
    /// Peptides x samples, a cell >= 0.5 meaning detected - the same convention
    /// <see cref="DetectionTest"/> uses.
    /// </param>
    public static IReadOnlyList<DetectionPairedRow> Run(
        double[,] detectionMatrix,
        IReadOnlyList<string> peptideIds,
        IReadOnlyList<SamplePair> pairs,
        MultipleTesting correction = MultipleTesting.BenjaminiHochberg)
    {
        var nPeptides = detectionMatrix.GetLength(0);
        if (peptideIds.Count != nPeptides)
            throw new ArgumentException("peptideIds length must match the matrix rows", nameof(peptideIds));
        if (pairs.Count == 0)
            throw new ArgumentException("no matched pairs", nameof(pairs));

        var rows = new DetectionPairedRow[nPeptides];
        var pValues = new double[nPeptides];

        for (var pep = 0; pep < nPeptides; pep++)
        {
            int detA = 0, detB = 0, onlyA = 0, onlyB = 0;
            foreach (var pair in pairs)
            {
                var inA = detectionMatrix[pep, pair.AColumn] >= 0.5;
                var inB = detectionMatrix[pep, pair.BColumn] >= 0.5;
                if (inA)
                    detA++;
                if (inB)
                    detB++;
                if (inA && !inB)
                    onlyA++;
                else if (inB && !inA)
                    onlyB++;
            }

            var p = McNemar.TwoSidedP(onlyA, onlyB);
            pValues[pep] = p;
            rows[pep] = new DetectionPairedRow(
                peptideIds[pep], pairs.Count, detA, detB, onlyA, onlyB,
                detA / (double)pairs.Count, detB / (double)pairs.Count, p, double.NaN);
        }

        var q = Fdr.Adjust(pValues, correction);
        for (var pep = 0; pep < nPeptides; pep++)
            rows[pep] = rows[pep] with { Q = q[pep] };

        return rows.OrderBy(r => r.P).ToList();
    }
}
