using System;
using System.IO;
using System.Runtime.CompilerServices;
using SkylinePrism.Core.Qc;

namespace SkylinePrism.Tests.TestSupport;

/// <summary>
/// Process-wide settings every test needs, applied once when the test assembly loads.
/// </summary>
/// <remarks>
/// The CLI resolves marker panels from the user's saved lists (<see cref="ProteinListSet.DefaultPath"/>)
/// plus the shipped ones. Left alone, a test naming "EV markers (core)" would read whatever the
/// developer running it has saved under that name - a saved list of the same name shadows the shipped
/// panel - so the same test would compare different things on different machines. Pointing the lists
/// file at a path that does not exist leaves the shipped panels alone, everywhere.
/// </remarks>
internal static class HermeticEnvironment
{
    [ModuleInitializer]
    internal static void UseOnlyShippedProteinLists() =>
        Environment.SetEnvironmentVariable(ProteinListSet.ListsPathVariable,
            Path.Combine(Path.GetTempPath(), $"prism-tests-no-saved-lists-{Environment.ProcessId}.json"));
}
