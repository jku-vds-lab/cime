/// <reference types="react" />
import { Dataset, DatasetType, IVector, PSEPlugin } from 'projection-space-explorer';
export declare class ChemPlugin extends PSEPlugin {
    type: DatasetType;
    createFingerprint(dataset: Dataset, vectors: IVector[], scale: number, aggregate: boolean): JSX.Element;
}
//# sourceMappingURL=ChemPlugin.d.ts.map