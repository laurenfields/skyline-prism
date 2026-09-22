using System;
using System.Collections.Generic;
using System.Linq;

namespace SkylinePrism.Core.DifferentialAnalysis;

/// <summary>
/// A few names and a count of the rest, for a status line or a note that would otherwise carry
/// hundreds of them.
/// </summary>
/// <remarks>
/// How many to show is the caller's: a status line has room for three subject ids, while a
/// "not detected" note is meant to be scanned and earns a dozen. What is NOT the caller's is the
/// decision to truncate at all - joining an unbounded list into a wrapping text block pushes the
/// plot it belongs to off the pane. Ticking all 65 shipped marker panels produced a 1,273-name note
/// exactly that way.
/// </remarks>
public static class NamePreview
{
    /// <summary>
    /// The first <paramref name="show"/> names, then <c>+N more</c>. Shows all of them when that is
    /// all there are, so the common small case reads naturally and never says "+0 more".
    /// </summary>
    public static string Of(IReadOnlyList<string> names, int show)
    {
        if (show < 1)
            throw new ArgumentOutOfRangeException(nameof(show), show, "Show at least one name.");

        return names.Count <= show
            ? string.Join(", ", names)
            : string.Join(", ", names.Take(show)) + $", +{names.Count - show} more";
    }
}
