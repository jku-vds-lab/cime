import * as React from "react";
import { DatasetType, IVector, PSEPlugin } from "projection-space-explorer";
import { ChemLegendParent } from "./ChemDetail/ChemDetail";

export class ChemPlugin extends PSEPlugin {
  type = DatasetType.Chem;

  createFingerprint(
    vectors: IVector[],
    scale: number,
    aggregate: boolean
  ): JSX.Element {
    return (
      <ChemLegendParent
        selection={vectors}
        aggregate={aggregate}
      ></ChemLegendParent>
    );
  }
}
